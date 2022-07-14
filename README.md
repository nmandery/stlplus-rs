# stlplus-rs

[![Latest Version](https://img.shields.io/crates/v/stlplus-rs.svg)](https://crates.io/crates/stlplus-rs)
[![Documentation](https://docs.rs/stlplus-rs/badge.svg)](https://docs.rs/stlplus-rs)
![ci](https://github.com/nmandery/stlplus-rs/workflows/CI/badge.svg)
[![dependency status](https://deps.rs/repo/github/nmandery/stlplus-rs/status.svg)](https://deps.rs/repo/github/nmandery/stlplus-rs)

Port of the enhanced Seasonal Trend Decomposition using Loess (STL) implementation from
https://github.com/hafen/stlplus (and the java implementation https://github.com/ruananswer/twitter-anomalyDetection-java)
to rust.

At the current stage this project is mostly a rough port of the above repositories, there surely is room for improvement - PRs welcome.

## Example

See [crates_io_downloads.rs](examples/crates_io_downloads.rs).

```rust
cargo run --example crates_io_downloads
```

![](example.png)


## License

<sup>
Licensed under either of <a href="LICENSE-APACHE">Apache License, Version
2.0</a> or <a href="LICENSE-MIT">MIT license</a> at your option.
</sup>
