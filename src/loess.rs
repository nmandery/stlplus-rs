use crate::{Degree, Error, NextOddInt, ValidateNotZero};
use itertools::izip;
use num_traits::{AsPrimitive, Float};

struct CLoessOutput<T> {
    result: Vec<T>,
    slopes: Vec<T>,
}

#[allow(clippy::too_many_arguments)]
fn c_loess<T>(
    xx: &[i64],
    yy: &[T],
    degree: Degree,
    span: usize,
    ww: &[T],
    m: &[usize],
    l_idx: &[usize],
    max_dist: &[usize],
) -> Result<CLoessOutput<T>, Error>
where
    T: Float + 'static,
    usize: AsPrimitive<T>,
    i64: AsPrimitive<T>,
{
    let n = xx.len();
    let n_m = m.len();

    let mut x = vec![T::zero(); span];
    let mut w = vec![T::zero(); span];
    let mut xw = vec![T::zero(); span];
    let mut x2w = vec![T::zero(); span];
    let mut x3w = vec![T::zero(); span];

    let mut result = vec![T::zero(); n_m];
    let mut slopes = vec![T::zero(); n_m];

    //let span3 = span;
    let span = if span > n { n } else { span };
    //let span2 = span.saturating_sub(1) / 2;

    //let offset = m.first().copied().expect("m.len() must be > 0");

    // loop through all values of m
    for i in 0..n_m {
        let mut a = T::zero();

        // get weights, x, and a
        for j in 0..span {
            // w[j] = T::zero();
            //x[j] = <i64 as AsPrimitive<T>>::as_(xx[l_idx[i] + j]) - m[i].as_();
            x[j] = xx[l_idx[i] + j].as_() - m[i].as_();

            let r = x[j].abs();
            let tmp1 = r / max_dist[i].validate_not_zero()?.as_();
            let tmp2 = T::one() - tmp1.powi(3);
            w[j] = tmp2.powi(3);

            // scale by user-defined weights
            w[j] = w[j] * ww[l_idx[i] + j];

            a = a + w[j];
        }

        if degree == Degree::Degree0 {
            let a1 = T::one() / a.validate_not_zero()?;
            for j in 0..span {
                result[i] = result[i] + w[j] * a1 * yy[l_idx[i] + j];
            }
        } else {
            // get xw, x2w, b, c for degree 1 or 2
            let mut b = T::zero();
            let mut c = T::zero();
            for j in 0..span {
                xw[j] = x[j] * w[j];
                x2w[j] = x[j] * xw[j];
                b = b + xw[j];
                c = c + x2w[j];
            }

            if degree == Degree::Degree1 {
                let det = T::one() / (a * c - b * b).validate_not_zero()?;
                let a1 = c * det;
                let b1 = -b * det;
                let c1 = a * det;
                for j in 0..span {
                    result[i] = result[i] + (w[j] * a1 + xw[j] * b1) * yy[l_idx[i] + j];
                    slopes[i] = slopes[i] + (w[j] * b1 + xw[j] * c1) * yy[l_idx[i] + j];
                }
            } else if degree == Degree::Degree2 {
                // get x3w, d, and e for degree 2
                let mut d = T::zero();
                let mut e = T::zero();
                for j in 0..span {
                    x3w[j] = x[j] * x2w[j];
                    d = d + x3w[j];
                    e = e + x3w[j] * x[j];
                }
                let mut a1 = e * c - d * d;
                let mut b1 = c * d - e * b;
                let mut c1 = b * d - c * c;
                let mut a2 = c * d - e * b;
                let mut b2 = e * a - c * c;
                let mut c2 = b * c - d * a;

                let det = T::one() / (a * a1 + b * b1 + c * c1).validate_not_zero()?;
                a1 = a1 * det;
                b1 = b1 * det;
                c1 = c1 * det;
                a2 = a2 * det;
                b2 = b2 * det;
                c2 = c2 * det;
                for j in 0..span {
                    result[i] =
                        result[i] + (w[j] * a1 + xw[j] * b1 + x2w[j] * c1) * yy[l_idx[i] + j];
                    slopes[i] =
                        slopes[i] + (w[j] * a2 + xw[j] * b2 + x2w[j] * c2) * yy[l_idx[i] + j];
                }
            } else {
                unreachable!()
            }
        }
    }
    Ok(CLoessOutput { result, slopes })
}

/// `x` and `y` are expected to be of the same length
pub fn loess_stl<VALUE>(
    x: &[i64],
    y: &[VALUE],
    span: usize,
    degree: Degree,
    m: &[usize],
    weights: &[VALUE],
    jump: usize,
) -> Result<Vec<VALUE>, Error>
where
    VALUE: Float + 'static + Copy,
    usize: AsPrimitive<VALUE>,
    i64: AsPrimitive<VALUE>,
{
    let n = y.len();
    let span = span.next_odd();
    let s2 = (span + 1) / 2;

    let (l_idx, r_idx) = if n < span {
        let l_idx = vec![0usize; m.len()];
        let r_idx = vec![n, m.len()];
        (l_idx, r_idx)
    } else {
        let mut count_small_than_s2 = 0usize;
        let mut count_large_than_s2 = 0usize;
        let mut count_large_than_n_minus_s2 = 0usize;
        for m_value in m {
            if *m_value < s2 {
                count_small_than_s2 += 1;
            } else if *m_value <= (n - s2) {
                count_large_than_s2 += 1;
            } else {
                count_large_than_n_minus_s2 += 1;
            }
        }

        let mut l_idx = vec![0usize; m.len()];
        let mut r_idx = vec![span.saturating_sub(1); m.len()];

        for i in count_small_than_s2..(count_small_than_s2 + count_large_than_s2) {
            l_idx[i] = m[i] - s2;
            r_idx[i] += l_idx[i];
        }
        for i in (count_small_than_s2 + count_large_than_s2)
            ..(count_small_than_s2 + count_large_than_s2 + count_large_than_n_minus_s2)
        {
            l_idx[i] = n.saturating_sub(span);
            r_idx[i] += l_idx[i];
        }
        (l_idx, r_idx)
    };

    let mut max_dist = vec![0usize; m.len()];
    izip!(max_dist.iter_mut(), m.iter().enumerate()).for_each(|(max_dist_value, (i, m_value))| {
        let aa = (*m_value as i64 - x[l_idx[i]] as i64).unsigned_abs() as usize;
        let bb = (x[r_idx[i]] as i64 - *m_value as i64) as usize;
        *max_dist_value = aa.max(bb);
    });

    if span > n {
        let add_value = (span - n) / 2;
        max_dist
            .iter_mut()
            .for_each(|md_value| *md_value += add_value);
    }
    let out = c_loess(x, y, degree, span, weights, m, &l_idx, &max_dist)?;

    // do interpolation
    if jump > 1 {
        //let at: Vec<VALUE> = (1..=n).map(|v| v.as_()).collect();
        let at: Vec<VALUE> = (0..n).map(|v| v.as_()).collect();
        Ok(interp(m, &out.result, &out.slopes, at.as_slice()))
    } else {
        Ok(out.result)
    }
}

/// This function was called `c_interp` in stlplus
fn interp<T>(m: &[usize], fits: &[T], slopes: &[T], at: &[T]) -> Vec<T>
where
    T: Float + 'static,
    usize: AsPrimitive<T>,
{
    let mut ans = vec![T::zero(); at.len()];

    let mut j = 0;
    izip!(at.iter(), ans.iter_mut()).for_each(|(at_value, ans_value)| {
        if *at_value > m[j + 1].as_() {
            j += 1;
        }

        let h: T = (m[j + 1] - m[j]).as_();
        let u = (*at_value - m[j].as_()) / h;
        let u2 = u.powi(2);
        let u3 = u2 * u;
        *ans_value = (2usize.as_() * u3 - 3usize.as_() * u2 + 1usize.as_()) * fits[j]
            + (3usize.as_() * u2 - 2usize.as_() * u3) * fits[j + 1]
            + (u3 - 2usize.as_() * u2 + u) * slopes[j] * h
            + (u3 - u2) * slopes[j + 1] * h;
    });
    ans
}
