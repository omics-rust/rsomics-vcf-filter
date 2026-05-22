use std::io::{self, Read};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Default)]
pub struct FilterConfig {
    pub min_qual: Option<f32>,
    pub pass_only: bool,
}

pub struct FilterStats {
    pub total: u64,
    pub passed: u64,
}

/// Filter VCF records by FILTER status and/or QUAL, passing kept lines through
/// verbatim. A tab byte-scan reads only the QUAL (col 5) and FILTER (col 6)
/// columns — no full record parse or re-serialization — so kept records keep
/// their exact original bytes (header included).
pub fn filter_vcf(
    input: &Path,
    output: &mut dyn io::Write,
    cfg: &FilterConfig,
) -> Result<FilterStats> {
    let raw = std::fs::read(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let data = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut d = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut d)
            .map_err(RsomicsError::Io)?;
        d
    } else {
        raw
    };

    let mut stats = FilterStats {
        total: 0,
        passed: 0,
    };
    for raw_line in data.split(|&b| b == b'\n') {
        let line = match raw_line.last() {
            Some(b'\r') => &raw_line[..raw_line.len() - 1],
            _ => raw_line,
        };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' {
            output.write_all(line).map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }
        stats.total += 1;

        let mut cols = line.split(|&b| b == b'\t');
        let qual = cols.nth(5);
        let filter = cols.next();

        if cfg.pass_only {
            let f = filter.unwrap_or(b".");
            if f != b"PASS" && f != b"." {
                continue;
            }
        }
        if let Some(min_q) = cfg.min_qual {
            // A present QUAL below the threshold is dropped; missing QUAL (".")
            // carries no score and is kept (matches noodles-era behaviour).
            if let Some(q) = qual.filter(|&q| q != b".") {
                let val: f32 = std::str::from_utf8(q)
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| {
                        RsomicsError::InvalidInput(format!(
                            "bad QUAL {:?}",
                            String::from_utf8_lossy(q)
                        ))
                    })?;
                if val < min_q {
                    continue;
                }
            }
        }

        output.write_all(line).map_err(RsomicsError::Io)?;
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
        stats.passed += 1;
    }

    Ok(stats)
}
