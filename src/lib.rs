use std::io::{self, Read};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_vcf_expr::{EvalContext, Expr};

#[derive(Default)]
pub struct FilterConfig {
    pub min_qual: Option<f32>,
    pub pass_only: bool,
    /// Include expression: keep records where expr is true (from `-i`).
    pub include_expr: Option<Expr>,
    /// Exclude expression: keep records where expr is false (from `-e`).
    pub exclude_expr: Option<Expr>,
}

pub struct FilterStats {
    pub total: u64,
    pub passed: u64,
}

/// Filter VCF records by FILTER status, QUAL, and optional bcftools-style
/// include/exclude expressions.
///
/// A tab byte-scan reads only the QUAL (col 5) and FILTER (col 6) columns
/// for the basic checks — no full record parse for kept records — so their
/// original bytes are preserved verbatim.  The expression engine reads the
/// full line when `-i/-e` is in use.
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

    let eval_ctx: Option<EvalContext> = match (&cfg.include_expr, &cfg.exclude_expr) {
        (Some(expr), None) => Some(EvalContext::new(expr.clone(), false)),
        (None, Some(expr)) => Some(EvalContext::new(expr.clone(), true)),
        (None, None) => None,
        (Some(_), Some(_)) => {
            return Err(RsomicsError::InvalidInput(
                "only one of -i or -e may be specified".into(),
            ));
        }
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
            // bcftools `-f PASS,.` semantics: keep the record when ANY of the
            // semicolon-separated FILTER tags is PASS or the missing-value dot.
            let keep = f
                .split(|&b| b == b';')
                .any(|tag| tag == b"PASS" || tag == b".");
            if !keep {
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
        if let Some(ctx) = &eval_ctx {
            let line_str = std::str::from_utf8(line)
                .map_err(|_| RsomicsError::InvalidInput("non-UTF-8 VCF line".into()))?;
            let result = ctx
                .eval_line(line_str, 0)
                .map_err(|e| RsomicsError::InvalidInput(format!("expression eval error: {e}")))?;
            // For site-level filtering: if any sample passes, keep the record.
            // For INFO/QUAL expressions the result is uniform across samples.
            // A record with no samples: result.pass has one element.
            let site_pass = result.pass.iter().any(|&p| p);
            if !site_pass {
                continue;
            }
        }

        output.write_all(line).map_err(RsomicsError::Io)?;
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
        stats.passed += 1;
    }

    Ok(stats)
}
