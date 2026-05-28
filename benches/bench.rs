use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vcf_filter(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-filter");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vcf = manifest.join("tests/golden/two.vcf");
    c.bench_function("rsomics-vcf-filter golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([vcf.to_str().unwrap(), "-o", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_filter);
criterion_main!(benches);
