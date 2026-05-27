use std::path::PathBuf;
use std::process::{Command, Stdio};
fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-filter"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bcftools_available() -> bool {
    Command::new("bcftools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Sorted CHROM\tPOS of data (non-header) records — identifies the kept set
/// independent of header/record reformatting differences between tools.
fn kept_positions(vcf: &[u8]) -> Vec<(String, String)> {
    let mut v: Vec<(String, String)> = String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| {
            let mut c = l.split('\t');
            Some((c.next()?.to_owned(), c.next()?.to_owned()))
        })
        .collect();
    v.sort();
    v
}

// --pass-only must keep the same record set as `bcftools view -f PASS,.`.
#[test]
fn pass_only_matches_bcftools() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }
    let vcf = fixture("two.vcf");
    let ours_out = Command::new(ours())
        .arg("--pass-only")
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(ours_out.status.success());
    let theirs = Command::new("bcftools")
        .args(["view", "-f", "PASS,."])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(theirs.status.success());
    assert_eq!(
        kept_positions(&ours_out.stdout),
        kept_positions(&theirs.stdout)
    );
}
// FILTER fields with semicolon-separated tags (e.g. "PASS;LowQual") must be kept
// by --pass-only when ANY tag is PASS — matching `bcftools view -f PASS,.`.
// The original implementation checked the whole field as a string literal, so
// "PASS;LowQual" != "PASS" and was wrongly dropped.
#[test]
fn pass_only_compound_filter_matches_bcftools() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }
    let vcf = fixture("mixed_filters.vcf");
    let ours_out = Command::new(ours())
        .arg("--pass-only")
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(ours_out.status.success());
    let theirs = Command::new("bcftools")
        .args(["view", "-f", "PASS,."])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(theirs.status.success());
    assert_eq!(
        kept_positions(&ours_out.stdout),
        kept_positions(&theirs.stdout),
        "compound FILTER like PASS;LowQual must be kept when PASS is one of the tags"
    );
}

#[test]
fn runs() {
    let out = Command::new(ours())
        .arg(fixture("two.vcf"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
