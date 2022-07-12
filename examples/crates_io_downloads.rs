use num_traits::Float;
use plotters::prelude::*;
use stlplus_rs::{stl_decompose, STLOptions, STLOutput};

#[rustfmt::skip]
fn daily_downloads() -> Vec<f32> {
    // see fetch_crates_io_testdata.py
    vec![
        108113.0, 103580.0, 85339.0, 40711.0, 45263.0, 84316.0, 108669.0, 110202.0, 106642.0,
        105171.0, 43743.0, 45082.0, 110190.0, 110791.0, 107311.0, 109093.0, 99330.0, 43257.0,
        47701.0, 113499.0, 111596.0, 104487.0, 103415.0, 109459.0, 47813.0, 42045.0, 118689.0,
        114456.0, 112742.0, 108783.0, 89554.0, 39153.0, 41106.0, 101895.0, 102104.0, 105788.0,
        //f32::NAN,
        110034.0, 107436.0, 47510.0, 43332.0, 112624.0, 116730.0, 109939.0, 94388.0, 90817.0,
        41149.0, 39748.0, 80703.0, 103159.0, 115283.0, 102565.0, 101867.0, 42052.0, 41680.0,
        111712.0, 119406.0, 112199.0, 120700.0, 109765.0, 46288.0, 41949.0, 126540.0, 128660.0,
        119960.0, 116728.0, 84615.0, 41455.0, 32614.0, 86439.0, 95445.0, 92745.0, 92750.0, 84713.0,
        38181.0, 35594.0, 92890.0, 98724.0, 88924.0, 87921.0, 91981.0, 30671.0, 29672.0, 70103.0,
        83425.0, 89964.0, 88415.0, 76165.0, 30552.0, 25771.0, 43555.0,
    ]
}

const OUT_FILE_NAME: &str = "downloads.png";

fn main() {
    let data = daily_downloads();

    let options = STLOptions {
        num_obs_per_period: 7,
        ..Default::default()
    };
    let decomp = stl_decompose(&data, &options).unwrap();

    plot(&data, &decomp, OUT_FILE_NAME);
}

fn plot<A>(input: &[f32], decomp: &STLOutput<f32>, out_file: A)
where
    A: AsRef<std::path::Path>,
{
    let root = BitMapBackend::new(OUT_FILE_NAME, (1000, 1000)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let root = root.margin(10, 10, 10, 10);
    let drawing_areas = root.split_evenly((4, 1));

    let (_min_value, max_value) = min_max(input);
    let (min_value, _max_value) = min_max(&decomp.seasonal);

    let sub_plots = vec![
        ("Input data", &drawing_areas[0], &BLACK, input),
        ("Seasonality", &drawing_areas[1], &BLUE, &decomp.seasonal),
        ("Trend", &drawing_areas[2], &BLUE, &decomp.trend),
        ("Remainder", &drawing_areas[3], &BLUE, &decomp.remainder),
    ];

    for (title, drawing_area, color, data_slice) in sub_plots.iter() {
        let mut chart = ChartBuilder::on(drawing_area)
            .set_label_area_size(LabelAreaPosition::Left, 80)
            .set_label_area_size(LabelAreaPosition::Bottom, 30)
            .build_cartesian_2d(0.0..(data_slice.len() as f32), min_value..max_value)
            .unwrap();

        chart
            .configure_mesh()
            .y_desc(*title)
            .x_desc("Time step")
            .draw()
            .unwrap();

        // plot all segments without a nan value
        let mut i = 0;
        let mut j = 0;
        while j + 1 < data_slice.len() {
            while i < data_slice.len() && data_slice[i].is_nan() {
                i += 1;
            }
            j = i;
            while j < data_slice.len() && !data_slice[j].is_nan() {
                j += 1;
            }
            if i != j {
                chart
                    .draw_series(LineSeries::new(
                        data_slice[i..j]
                            .iter()
                            .enumerate()
                            .map(|(idx, value)| ((idx + i) as f32, *value)),
                        color,
                    ))
                    .unwrap();
            }
            i = j;
        }
    }

    root.present().unwrap();
    println!(
        "Plot has been saved to {}",
        out_file.as_ref().to_str().unwrap()
    );
}

fn min_max<T>(slice: &[T]) -> (T, T)
where
    T: Float,
{
    let (min_value, max_value) = slice.iter().fold(
        (*slice.first().unwrap(), *slice.first().unwrap()),
        |acc, value| {
            (
                if *value < acc.0 { *value } else { acc.0 },
                if *value > acc.1 { *value } else { acc.1 },
            )
        },
    );
    (min_value, max_value)
}
