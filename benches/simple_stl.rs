use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use stlplus_rs::{stl_decompose, STLOptions};

#[rustfmt::skip]
fn daily_downloads() -> Vec<f32> {
    // see fetch_crates_io_testdata.py
    vec![
        108113.0, 103580.0, 85339.0, 40711.0, 45263.0, 84316.0, 108669.0, 110202.0, 106642.0,
        105171.0, 43743.0, 45082.0, 110190.0, 110791.0, 107311.0, 109093.0, 99330.0, 43257.0,
        47701.0, 113499.0, 111596.0, 104487.0, 103415.0, 109459.0, 47813.0, 42045.0, 118689.0,
        114456.0, 112742.0, 108783.0, 89554.0, 39153.0, 41106.0, 101895.0, 102104.0, 105788.0,
        110034.0, 107436.0, 47510.0, 43332.0, 112624.0, 116730.0, 109939.0, 94388.0, 90817.0,
        41149.0, 39748.0, 80703.0, 103159.0, 115283.0, 102565.0, 101867.0, 42052.0, 41680.0,
        111712.0, 119406.0, 112199.0, 120700.0, 109765.0, 46288.0, 41949.0, 126540.0, 128660.0,
        119960.0, 116728.0, 84615.0, 41455.0, 32614.0, 86439.0, 95445.0, 92745.0, 92750.0, 84713.0,
        38181.0, 35594.0, 92890.0, 98724.0, 88924.0, 87921.0, 91981.0, 30671.0, 29672.0, 70103.0,
        83425.0, 89964.0, 88415.0, 76165.0, 30552.0, 25771.0, 43555.0,
    ]
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("simple_stl");
    group.sample_size(1000);
    group.warm_up_time(Duration::from_secs(2));
    group.bench_function("simple STL", |bencher| {
        let data = daily_downloads();

        let options = STLOptions {
            num_obs_per_period: 7,
            ..Default::default()
        };
        bencher.iter(|| {
            let _decomp = stl_decompose(&data, &options).unwrap();
        });
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
