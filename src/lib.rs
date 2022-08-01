#[cfg(feature = "arrow2")]
pub mod arrow;
pub mod error;
mod loess;
mod util;

use crate::error::Error;
use crate::loess::loess_stl;
use crate::util::{NextOddInt, SumAgg, ValidateNotZero};
use itertools::izip;
use num_traits::{AsPrimitive, Float};
use std::fmt::Debug;
pub use util::NextOdd;

#[derive(Debug, Eq, PartialEq, Copy, Clone)]
pub enum Degree {
    Degree0,
    Degree1,
    Degree2,
}

impl Degree {
    /// index in the COEF_XX static arrays
    fn coef_index(&self) -> usize {
        match self {
            Degree::Degree0 => 0,
            Degree::Degree1 => 1,
            Degree::Degree2 => 2,
        }
    }
}

macro_rules! impl_tryfrom_int_for_degree {
    ($($t:ty),*) => {
        $(
            impl TryFrom<$t> for Degree {
                type Error = Error;

                fn try_from(value: $t) -> Result<Self, Self::Error> {
                    match value {
                        1 => Ok(Self::Degree0),
                        2 => Ok(Self::Degree1),
                        3 => Ok(Self::Degree2),
                        _ => Err(Error::InvalidDegree),
                    }
                }
            }
        )*
    };
}
impl_tryfrom_int_for_degree!(i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize);

///
/// Comments are still from the original R implementation and Java port.
pub struct STLOptions {
    /// The number of observations in each cycle of the seasonal component, n_p
    pub num_obs_per_period: usize,

    /// s.window either the character string \code{"periodic"} or the span (in lags) of the loess window for seasonal extraction,
    /// which should be odd.  This has no default.
    ///
    /// None is used for periodic
    pub s_window: Option<usize>,

    /// s.degree degree of locally-fitted polynomial in seasonal extraction.  Should be 0, 1, or 2.
    pub s_degree: Degree,

    /// t.window the span (in lags) of the loess window for trend extraction, which should be odd.
    /// If \code{NULL}, the default, \code{nextodd(ceiling((1.5*period) / (1-(1.5/s.window))))}, is taken.
    pub t_window: Option<usize>,

    /// t.degree degree of locally-fitted polynomial in trend extraction.  Should be 0, 1, or 2.
    pub t_degree: Degree,

    /// l.window the span (in lags) of the loess window of the low-pass filter used for each subseries.
    ///
    /// Defaults to the smallest odd integer greater than or equal to \code{n.p}
    /// which is recommended since it prevents competition between the trend and seasonal components.
    /// If not an odd integer its given value is increased to the next odd one.
    pub l_window: Option<usize>,

    /// l.degree degree of locally-fitted polynomial for the subseries low-pass filter.  Should be 0, 1, or 2.
    pub l_degree: Degree,

    /// s.jump s.jump,t.jump,l.jump,fc.jump integers at least one to increase speed of the respective smoother.
    /// Linear interpolation happens between every \code{*.jump}th value.
    pub s_jump: Option<usize>,

    /// t.jump
    pub t_jump: Option<usize>,

    /// l.jump
    pub l_jump: Option<usize>,

    /// critfreq the critical frequency to use for automatic calculation of smoothing windows for the trend and high-pass filter.
    pub critfreq: f64,

    /// The number of passes through the inner loop, n_i
    pub number_of_inner_loop_passes: u32,

    /// The number of robustness iterations of the outer loop, n_o
    pub number_of_robustness_iterations: u32,
}

impl Default for STLOptions {
    fn default() -> Self {
        let num_obs_per_period = 4;
        Self {
            num_obs_per_period,
            s_window: None,
            s_degree: Degree::Degree1,
            t_window: None,
            t_degree: Degree::Degree1,
            l_window: None,
            l_degree: Degree::Degree1,
            s_jump: None,
            t_jump: None,
            l_jump: None,
            critfreq: 0.05,
            number_of_inner_loop_passes: 2,
            number_of_robustness_iterations: 1,
            //number_of_data_points: num_obs_per_period * 2,
        }
    }
}

pub struct STLOutput<VALUE> {
    pub trend: Vec<VALUE>,
    pub seasonal: Vec<VALUE>,
    pub remainder: Vec<VALUE>,
}

pub fn stl_decompose<VALUE>(
    values: &[VALUE],
    options: &STLOptions,
) -> Result<STLOutput<VALUE>, Error>
where
    VALUE: Float + 'static + Copy,
    usize: AsPrimitive<VALUE>,
    i64: AsPrimitive<VALUE>,
{
    let n = values.len();
    let times_i64: Vec<_> = (0..n).map(|v| v as i64).collect();

    stl_decompose_with_time(&times_i64, values, options)
}

pub fn stl_decompose_with_time<TIME, VALUE>(
    times: &[TIME],
    values: &[VALUE],
    options: &STLOptions,
) -> Result<STLOutput<VALUE>, Error>
where
    VALUE: Float + 'static + Copy,
    usize: AsPrimitive<VALUE>,
    TIME: AsPrimitive<i64>,
    i64: AsPrimitive<VALUE>,
{
    if values.is_empty() || times.is_empty() {
        return Err(Error::EmptyInputSlice);
    }
    let n = values.len();
    if times.len() != n {
        return Err(Error::InputSlicesDifferingLength);
    }
    let times_i64: Vec<i64> = times.iter().map(|t| t.as_()).collect();

    let validated_options = validate_options(options, n)?;
    let mut trend = vec![VALUE::zero(); n];
    let mut seasonal = vec![VALUE::zero(); n];
    let mut deseasonalized = vec![VALUE::zero(); n];

    // cycleSubIndices will keep track of what part of the seasonal each observation belongs to
    let cycle_sub_indices: Vec<_> = (1..=validated_options.num_obs_per_period)
        .cycle()
        .take(values.len())
        .collect();
    let weights = vec![VALUE::one(); n];
    let mut detrend = vec![VALUE::zero(); n];
    // todo: missing stuff from java impl?

    let mut cycle_sub = Vec::with_capacity(
        (n as f64 / validated_options.num_obs_per_period as f64).ceil() as usize / 2,
    );
    let mut sub_weights = Vec::with_capacity(cycle_sub.capacity());

    let (cs1, cs2) = {
        let mut cs1 = Vec::with_capacity(validated_options.num_obs_per_period);
        let mut cs2 = Vec::with_capacity(validated_options.num_obs_per_period);
        for i in 0..validated_options.num_obs_per_period {
            cs1.push(cycle_sub_indices[i]);
            cs2.push(cycle_sub_indices[values.len() - validated_options.num_obs_per_period + i]);
        }
        (cs1, cs2)
    };

    let l_ev = Ev::new(n, validated_options.l_jump);
    let t_ev = Ev::new(n, validated_options.t_jump);

    let mut c = vec![VALUE::nan(); n + 2 * validated_options.num_obs_per_period];

    // start and end indices for after adding in extra n.p before and after
    let c_start_idx = validated_options.num_obs_per_period;
    let c_end_idx = n - 1 + validated_options.num_obs_per_period;

    for _outer_iteration_i in 1..=options.number_of_robustness_iterations {
        for _inner_iteration_i in 1..=options.number_of_inner_loop_passes {
            // Step 1: detrending
            izip!(detrend.iter_mut(), values.iter(), trend.iter()).for_each(|(dt, v, t)| {
                *dt = *v - *t;
            });

            // step 2: smoothing of cycle-subseries
            for i in 0..validated_options.num_obs_per_period {
                cycle_sub.clear();
                sub_weights.clear();
                let mut j = i;
                while j < n {
                    if cycle_sub_indices[j] == i + 1 {
                        cycle_sub.push(detrend[j]);
                        sub_weights.push(weights[j])
                    }
                    j += validated_options.num_obs_per_period;
                }

                let weight_mean_ans = weight_mean(&cycle_sub, &sub_weights)?;
                j = i;
                while j < validated_options.num_obs_per_period {
                    if cs1[j] == i + 1 {
                        c[j] = weight_mean_ans;
                    }
                    j += validated_options.num_obs_per_period;
                }
                j = i;
                while j < n {
                    if cycle_sub_indices[j] == i + 1 {
                        c[j + validated_options.num_obs_per_period] = weight_mean_ans;
                    }
                    j += validated_options.num_obs_per_period;
                }
                for j in 0..validated_options.num_obs_per_period {
                    if cs2[j] == i + 1 {
                        c[j + validated_options.num_obs_per_period + n] = weight_mean_ans;
                    }
                }
            }

            // Step 3: Low-pass filtering of collection of all the cycle-subseries
            // moving averages
            let ma3 = cycle_subseries_moving_averages(&c, validated_options.num_obs_per_period);

            // Step 4: Detrend smoothed cycle-subseries
            let l = loess_stl(
                &times_i64,
                &ma3,
                validated_options.l_window,
                validated_options.l_degree,
                l_ev.as_slice(),
                &weights,
                validated_options.l_jump,
            )?;

            // Step 5: Deseasonalize
            izip!(
                seasonal.iter_mut(),
                (&c)[c_start_idx..=c_end_idx].iter(),
                l.iter(),
                values.iter(),
                deseasonalized.iter_mut()
            )
            .for_each(|(s, c, l, v, d)| {
                *s = *c - *l;
                *d = *v - *s;
            });

            // Step 6: Trend Smoothing
            trend = loess_stl(
                &times_i64,
                &deseasonalized,
                validated_options.t_window,
                validated_options.t_degree,
                t_ev.as_slice(),
                &weights,
                validated_options.t_jump,
            )?;
        }
    }

    let remainder: Vec<_> = izip!(values.iter(), trend.iter(), seasonal.iter())
        .map(|(v, t, s)| *v - *t - *s)
        .collect();

    Ok(STLOutput {
        trend,
        seasonal,
        remainder,
    })
}

struct Ev {
    n: usize,
    array_min_len: usize,
    storage_vec: Vec<usize>,
}

impl Ev {
    fn new(n: usize, jump: usize) -> Self {
        let array_min_len = (n as f64 / jump as f64).ceil() as usize;
        let mut storage_vec = vec![0usize; array_min_len + 1];

        let mut i = 0;
        let mut j = 0;
        while i < array_min_len {
            storage_vec[i] = j + 1;
            i += 1;
            j += jump
        }

        // always have the last element == n
        storage_vec[array_min_len] = n;

        Self {
            n,
            array_min_len,
            storage_vec,
        }
    }

    /// return a slice where the last element == `n`
    fn as_slice(&self) -> &[usize] {
        if self.storage_vec[self.array_min_len - 1] != self.n {
            &self.storage_vec
        } else {
            &self.storage_vec[0..self.array_min_len]
        }
    }
}

struct ValidatedOptions {
    num_obs_per_period: usize,
    //number_of_data_points: usize,
    //s_window: usize,
    //s_degree: Degree,
    //s_jump: usize,
    t_window: usize,
    t_degree: Degree,
    t_jump: usize,
    l_window: usize,
    l_degree: Degree,
    l_jump: usize,
    //periodic: bool,
}

fn validate_options(options: &STLOptions, num_values: usize) -> Result<ValidatedOptions, Error> {
    let num_obs_per_period = if options.num_obs_per_period >= 4 {
        options.num_obs_per_period
    } else {
        return Err(Error::InvalidNumObsPerPeriod);
    };

    /*
    let number_of_data_points = if options.number_of_data_points > 2 * num_obs_per_period {
        options.number_of_data_points
    } else {
        return Err(Error::InvalidNumberOfDataPoints);
    };

     */

    let l_degree = options.l_degree;
    let l_window = options.l_window.unwrap_or(num_obs_per_period).next_odd();
    let l_jump = options.l_jump.unwrap_or_else(|| window_to_jump(l_window));

    let (s_window, s_degree, _s_jump, _periodic) = if let Some(s_window) = options.s_window {
        let s_window = validate_window(s_window)?;
        let s_jump = options.s_jump.unwrap_or_else(|| window_to_jump(s_window));
        (s_window, options.s_degree, s_jump, false)
    } else {
        // periodic
        let s_window = 10 * num_values + 1;
        let s_degree = Degree::Degree0;
        let s_jump = window_to_jump(s_window);
        (s_window, s_degree, s_jump, true)
    };

    let t_degree = options.t_degree;
    let t_window = if let Some(t_window) = options.t_window {
        validate_window(t_window)?
    } else {
        get_t_window(
            t_degree,
            s_degree,
            s_window,
            num_obs_per_period,
            options.critfreq,
        )?
    };
    let t_jump = options.t_jump.unwrap_or_else(|| window_to_jump(t_window));

    Ok(ValidatedOptions {
        num_obs_per_period,
        //number_of_data_points,
        //s_window,
        //s_degree,
        //s_jump,
        t_window,
        t_degree,
        t_jump,
        l_window,
        l_degree,
        l_jump,
        //periodic,
    })
}

static COEFS_A: [[f64; 2]; 2] = [
    [0.000103350651767650, 3.81086166990428e-6],
    [-0.000216653946625270, 0.000708495976681902],
];
static COEFS_B: [[f64; 2]; 3] = [
    [1.42686036792937, 2.24089552678906],
    [-3.1503819836694, -3.30435316073732],
    [5.07481807116087, 5.08099438760489],
];
static COEFS_C: [[f64; 2]; 3] = [
    [1.66534145060448, 2.33114333880815],
    [-3.87719398039131, -1.8314816166323],
    [6.46952900183769, 1.85431548427732],
];

fn get_t_window(
    t_degree: Degree,
    s_degree: Degree,
    s_window: usize,
    num_obs_per_period: usize,
    critfreq: f64,
) -> Result<usize, Error> {
    let s_index = s_degree.coef_index();
    let t_index = t_degree.coef_index();

    // estimate critical frequency for seasonal
    let betac0 = COEFS_A[1][s_index].mul_add(critfreq, COEFS_A[0][s_index]);
    let betac1 = COEFS_B[1][s_index].mul_add(critfreq, COEFS_B[0][s_index])
        + COEFS_B[2][s_index] * critfreq.powi(2);
    let betac2 = COEFS_C[1][s_index].mul_add(critfreq, COEFS_C[0][s_index])
        + COEFS_C[2][s_index] * critfreq.powi(2);

    let f_c = (1.0 - (betac0 + betac1 / s_window as f64 + betac2 / s_window.pow(2) as f64))
        / num_obs_per_period as f64;

    // choose
    let betat0 = COEFS_A[1][t_index].mul_add(critfreq, COEFS_A[0][t_index]);
    let betat1 = COEFS_B[1][t_index].mul_add(critfreq, COEFS_B[0][t_index])
        + COEFS_B[2][t_index] * critfreq.powi(2);
    let betat2 = COEFS_C[1][t_index].mul_add(critfreq, COEFS_C[0][t_index])
        + COEFS_C[2][t_index] * critfreq.powi(2);

    let betat00 = betat0 - f_c;

    Ok(
        ((-betat1 - (betat1.powi(2) - 4.0 * betat00 * betat2).sqrt()) / (2.0 * betat00)).next_odd()
            as usize,
    )
}

fn validate_window(window: usize) -> Result<usize, Error> {
    if window < 1 {
        Err(Error::InvalidWindow)
    } else {
        Ok(window)
    }
}

fn window_to_jump(window: usize) -> usize {
    (window as f64 / 10.0).ceil() as usize
}

/// cycle-subseries moving averages
///
/// `num_obs_per_period` is the periodicity `n_p`
///
/// This function was called `c_ma` in stlplus
pub(crate) fn cycle_subseries_moving_averages<F>(x: &[F], num_obs_per_period: usize) -> Vec<F>
where
    F: Float + 'static + Copy,
    usize: AsPrimitive<F>,
{
    let nn = x.len().saturating_sub(num_obs_per_period * 2);

    let mut ans = vec![F::zero(); x.len() - 2 * num_obs_per_period];
    let mut ma = vec![F::zero(); nn + num_obs_per_period + 1];
    let mut ma2 = vec![F::zero(); nn + 2];

    let mut ma_tmp = x[0..num_obs_per_period].sum_agg();
    ma[0] = ma_tmp / num_obs_per_period.as_();

    for i in num_obs_per_period..(nn + 2 * num_obs_per_period) {
        ma_tmp = ma_tmp - x[i - num_obs_per_period] + x[i];
        ma[i - num_obs_per_period + 1] = ma_tmp / num_obs_per_period.as_();
    }

    ma_tmp = (&ma[0..num_obs_per_period]).sum_agg();
    ma2[0] = ma_tmp / num_obs_per_period.as_();
    for i in num_obs_per_period..(nn + num_obs_per_period + 1) {
        ma_tmp = ma_tmp - ma[i - num_obs_per_period] + ma[i];
        ma2[i - num_obs_per_period + 1] = ma_tmp / num_obs_per_period.as_();
    }

    ma_tmp = (&ma2[0..3]).sum_agg();
    ans[0] = ma_tmp / 3usize.as_();
    for i in 3..(nn + 2) {
        ma_tmp = ma_tmp - ma2[i - 3] + ma2[i];
        ans[i - 2] = ma_tmp / 3usize.as_();
    }

    ans
}

fn weight_mean<T>(x: &[T], w: &[T]) -> Result<T, Error>
where
    T: Float,
{
    let (sum, sum_w) = x.iter().zip(w.iter()).fold(
        (T::zero(), T::zero()),
        |(sum, sum_w), (x_value, w_value)| {
            if !x_value.is_nan() {
                (sum + (*x_value * *w_value), sum_w + *w_value)
            } else {
                (sum, sum_w)
            }
        },
    );
    Ok(sum / sum_w.validate_not_zero()?)
}

#[cfg(test)]
mod tests {
    use crate::{stl_decompose, STLOptions};

    #[test]
    fn c_ma() {
        let input = vec![1.0f32, 1.0, 2.0, 2.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0];
        let n_p = 3;
        let out = super::cycle_subseries_moving_averages(&input, n_p);
        dbg!(out);
    }

    /// https://github.com/nmandery/stlplus-rs/issues/1
    #[test]
    fn indexing_within_bounds() {
        let input = vec![0.0f32; 2581];
        let options = STLOptions {
            num_obs_per_period: 365,
            ..Default::default()
        };
        // should not cause out-of-bounds panic
        stl_decompose(&input, &options).unwrap();
    }
}
