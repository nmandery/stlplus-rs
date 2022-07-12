use crate::{Error, STLOptions, STLOutput};
use arrow2::array::{Array, PrimitiveArray};
use arrow2::types::NativeType;
use num_traits::{AsPrimitive, Float};

pub struct ArrowDecomposeOutput<VALUE: NativeType> {
    pub trend: PrimitiveArray<VALUE>,
    pub seasonal: PrimitiveArray<VALUE>,
    pub remainder: PrimitiveArray<VALUE>,
}

impl<VALUE> From<STLOutput<VALUE>> for ArrowDecomposeOutput<VALUE>
where
    VALUE: NativeType + Float,
{
    fn from(decomp_output: STLOutput<VALUE>) -> Self {
        Self {
            trend: to_primitive_array(decomp_output.trend),
            seasonal: to_primitive_array(decomp_output.seasonal),
            remainder: to_primitive_array(decomp_output.remainder),
        }
    }
}

pub fn stl_decompose_arrow<VALUE>(
    values: PrimitiveArray<VALUE>,
    options: &STLOptions,
) -> Result<ArrowDecomposeOutput<VALUE>, Error>
where
    VALUE: Float + NativeType + 'static + Copy,
    usize: AsPrimitive<VALUE>,
    i64: AsPrimitive<VALUE>,
{
    let values_vec: Vec<_> = (0..values.len()).map(|i| get_or_nan(&values, i)).collect();
    super::stl_decompose(&values_vec, options).map(|decomp_out| decomp_out.into())
}

pub fn stl_decompose_arrow_with_time<TIME, VALUE>(
    times: PrimitiveArray<TIME>,
    values: PrimitiveArray<VALUE>,
    options: &STLOptions,
) -> Result<ArrowDecomposeOutput<VALUE>, Error>
where
    VALUE: Float + NativeType + 'static + Copy,
    TIME: NativeType + AsPrimitive<i64>,
    usize: AsPrimitive<VALUE>,
    i64: AsPrimitive<VALUE>,
{
    // all times must be set, invalid values are not allowed
    if !all_values_valid(&times) {
        return Err(Error::InputSliceInvalidValues);
    }

    let values_vec: Vec<_> = (0..values.len()).map(|i| get_or_nan(&values, i)).collect();
    let times: Vec<_> = (0..times.len()).map(|i| times.value(i)).collect();
    super::stl_decompose_with_time(&times, &values_vec, options).map(|decomp_out| decomp_out.into())
}

#[inline(always)]
fn get_or_nan<T>(pa: &PrimitiveArray<T>, position: usize) -> T
where
    T: Float + NativeType,
{
    if pa.is_valid(position) {
        pa.value(position)
    } else {
        T::nan()
    }
}

fn to_primitive_array<T>(mut contents: Vec<T>) -> PrimitiveArray<T>
where
    T: Float + NativeType,
{
    PrimitiveArray::from_iter(
        contents
            .drain(..)
            .map(|v| if v.is_nan() { None } else { Some(v) }),
    )
}

fn all_values_valid<T>(array: &PrimitiveArray<T>) -> bool
where
    T: NativeType,
{
    array
        .validity()
        .map(|validity| validity.null_count())
        .unwrap_or(0)
        == 0
}

#[cfg(test)]
mod tests {
    use arrow2::array::PrimitiveArray;

    #[test]
    fn all_values_valid() {
        let pa1 = PrimitiveArray::<i32>::from([Some(1), Some(2), Some(3)]);
        assert!(super::all_values_valid(&pa1));

        let pa2 = PrimitiveArray::<i32>::from([Some(1), None, Some(3)]);
        assert!(!super::all_values_valid(&pa2));
    }
}
