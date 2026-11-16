# lrcpq 0.3.0

* **Defaults now match the LRCQ and LRCP papers.** Results from code that relied on the old defaults will change; pass the old values to reproduce them.
  * `lrcq()`, `run_lrcq()` and `lrcq_window()` use `method = "ols"` (was `"wls"`).
  * `lrcq()` and `run_lrcq()` return skewness-corrected p-values (`p_gt0`, `p_gt1`) by default with OLS (new argument `skew`, default `NULL` = on with OLS).
  * `run_lrcq()` leaves estimates unrectified (`rectify = "none"`; was `"C"`).
  * `run_lrcp()` screens tags with W > 0 at LRCQ FDR q < 0.05 (new arguments `prescreen = "fdr"` and `q_max`). The old screen on `w_raw >= threshold` stays available as `prescreen = "threshold"`.
  * Command line: `lrcpq lrcq` uses OLS with skewness correction (`--no-skew` turns it off; `--skew` is no longer needed). `lrcpq lrcp` uses the FDR prescreen (`--prescreen threshold --threshold 10` for the old behaviour; `--q` sets the level).
  * To reproduce 0.2 results: `run_lrcq(..., method = "wls", rectify = "C", skew = FALSE)` and `run_lrcp(..., prescreen = "threshold", threshold = 10)`.

# lrcpq 0.2.1

* Citation: `citation("lrcpq")` now gives the LRCQ and LRCP bioRxiv preprints (DOIs 10.64898/2026.11.02.763731 and 10.64898/2026.11.02.763742). The README no longer says the repository is private.
* Fix: `run_lrcq()` (and so `*.lrcq.tsv.gz` from the command line) no longer writes a duplicate `block` column.
* Fix: command-line runs append to `OUT.log` instead of overwriting it, so `lrcq` and `lrcp` runs with the same `--out` keep one log.
* `pair_sum = "twice"` replaces `"chapter"` (still accepted as an alias). Documentation no longer refers to thesis chapters or slides; `lrcp_distal()` now cites theory supplement S4.6.
* Getting-started vignette: text updated to match a fresh run.

# lrcpq 0.2.0

* The package moved to its own repository, unfated/LRCPQ.
* New: `munge_sumstats()` and `build_z_matrix()` take raw GWAS files (LDSC column synonyms; Z from BETA/SE, OR/SE or signed P), then filter and align them to the LD reference.
* New: an `ldsc.py`-style command-line tool (`lrcpq_cli()`, `inst/cli/lrcpq`) with `munge`, `lrcq` and `lrcp` commands.
* New: a bundled chr22 example (50 simulated GWAS, 1000G EUR LD), used by the getting-started vignette.
* New: skew-corrected LRCQ p-values (`p_skew()`, `lrcq(..., skew = TRUE)`).
* New: a genome-wide LRCP screen with simulated-null calibration (`lrcp_screen()`).
* New: `check_scale()` checks that n, h2 and M are on one scale.
* Fix: `p_gt1` now tests against `w_null1`, the expected tag estimate with no enrichment. It previously used the clump size.
* Change: LRCQ model SEs (and third cumulants) now use the background-filled rule-A plug-in by default (`se_rectify = "bgA"`), matching the LRCQ paper (theory supplement S3.10a). The old default `"A"` drops the polygenic background, which makes SEs 10-25% too small under dense enrichment. It is still available, along with `"eb"` and `"one"`. Estimates do not change.
* New: the prebuilt 1000G EUR HapMap3 PLINK reference is attached to the GitHub release.
* `panukb_select()` returns total `n` (`n_eff` is kept as an alias).
* `lrcp_gene()` reports bootstrap SDs, and `check_gcov()` warns when the intercept matrix is not positive semi-definite.

# lrcpq 0.1.0

* First version: simulation engine, LD utilities, LD score regression, and the LRCQ and LRCP estimators.
