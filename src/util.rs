use crate::Error;
use num_traits::{AsPrimitive, Float, Zero};
use std::ops::Add;

pub trait NextOddInt {
    fn next_odd(&self) -> Self;
}

macro_rules! impl_nextoddint {
    ($($t:ty),*) => {
        $(
            impl NextOddInt for $t {
                fn next_odd(&self) -> Self {
                    if self % 2 == 0 {
                        self + 1
                    } else {
                        *self
                    }
                }
            }
        )*
    };
}
impl_nextoddint!(i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize);

pub trait NextOdd {
    fn next_odd(&self) -> i64;
}

impl<T> NextOdd for T
where
    T: Float + 'static + Copy + AsPrimitive<i64>,
{
    fn next_odd(&self) -> i64 {
        let rounded: i64 = self.round().as_();
        rounded.next_odd()
    }
}

pub trait ValidateNotZero {
    fn validate_not_zero(self) -> Result<Self, Error>
    where
        Self: Sized;
}

impl<T> ValidateNotZero for T
where
    T: Zero,
{
    fn validate_not_zero(self) -> Result<Self, Error> {
        if self.is_zero() {
            Err(Error::UnexpectedZero)
        } else {
            Ok(self)
        }
    }
}

pub trait SumAgg {
    type Output;

    fn sum_agg(&self) -> Self::Output;
}

impl<T> SumAgg for [T]
where
    T: Zero + Copy + Add,
{
    type Output = T;

    fn sum_agg(&self) -> Self::Output {
        self.iter().fold(T::zero(), |a, b| a + *b)
    }
}
