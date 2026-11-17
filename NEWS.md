# lrcpq 0.3.0

* **Defaults now match the LRCQ and LRCP papers.** Results from code that relied on the old defaults will change; pass the old values to reproduce them.
  * `lrcq()`, `run_lrcq()` and `lrcq_window()` use `method = "ols"` (was `"wls"`).
  * `lrcq()` and `run_lrcq()` return skewness-corrected p-values (`p_gt0`, `p_gt1`) by default with OLS (new argument `skew`, default `NULL` = on with OLS).
  * `run_lrcq()` leaves estimates unrectified (`rectify = "none"`; was `"C"`).
  * `run_lrcp()` screens tags with W > 0 at LRCQ FDR q < 0.05 (new arguments `prescreen = "fdr"` and `q_max`). The old screen on `w_raw >= threshold` stays available as `prescreen = "threshold"`.
  * Command line: `lrcpq lrcq` uses OLS with skewness correction (`--no-skew` turns it off; `--skew` is no longer needed). `lrcpq lrcp` uses the FDR prescreen (`--prescreen threshold --threshold 10` for the old behaviour; `--q` sets the level).
  * To reproduce 0.2 results: `run_lrcq(..., method = "wls", rectify = "C", skew = FALSE)` and `run_lrcp(..., prescreen = "threshold", threshold = 10)`.
* **Score-type LRCQ tag tests** (`test = "score"`, the new default with OLS; CLI `--wald` for the old test). The cumulants of each tag's estimate are evaluated with the tag's own enrichment at its null value (0 for `p_gt0`, `w_null1` for `p_gt1`) and every other SNP at the background-filled plug-in. The Wald-type test capped z near 1.7 for the largest estimates. p-values use the normal reference up to score z = 1.5, the shifted gamma up to 2, and the new kurtosis-matched chi-square `p_kurt()` (Lee, Wu and Lin 2012) above. New columns `se_gt0`, `z_gt0`, `se_gt1`, `z_gt1` and `test`; `se` stays the Wald-type SE. Matches the LRCQ paper's re-analysis. To reproduce 0.2 p-values add `test = "wald"`.
* **Adaptive bootstrap in `lrcp_gene()` and `run_lrcp()`.** A pair whose bootstrap p-value is near the floor 1/(n_boot + 1) (fewer than `boot_min_exceed = 10` exceedances) gets further draws, in rounds of 10 times the draws so far, up to `n_boot_max = 1e5` (CLI `--n-boot-max`). Strong pairs no longer tie at the floor, so their q-values are informative. Other pairs keep their p-values. The draws behind each p-value are reported (`n_boot_used`; column `n_boot` in `run_lrcp()`). `n_boot_max = n_boot` gives the 0.2 behaviour.
* `lrcp_gene()` no longer builds the (n1 n2)-by-(n1 n2) covariance of the gene-pair matrix: SEs, the max-|T| simulation and the Liu-Satterthwaite `p_Q` now use its Kronecker factors. Results are unchanged; memory is now O(n1^2 + n2^2) without a class null.

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
