#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("empty input slice given")]
    EmptyInputSlice,

    #[error("input slice of differing length")]
    InputSlicesDifferingLength,

    #[error("input slice contains invalid values")]
    InputSliceInvalidValues,

    #[error("Invalid degree - must be 1, 2 or 3")]
    InvalidDegree,

    #[error("Invalid window - must be > 0")]
    InvalidWindow,

    #[error("length of sub_labels must match num_obs_per_period")]
    InvalidSubLabelsLength,

    #[error("sub_labels must be unique")]
    SubLabelsMustBeUnique,

    #[error("num_obs_per_period is invalid")]
    InvalidNumObsPerPeriod,
    // InvalidNumberOfDataPoints,
    #[error("unexpected 0 encountered")]
    UnexpectedZero,
}
