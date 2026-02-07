# Project Notes

## Completed: Exact Bias-Corrected ell_t (Paper Eq 13-14)

### Status: DONE — all 8 steps implemented, tests pass, simulation validated

### Context
The current code uses a raw approximation for the observation-level estimating equation:
- `ell_t = e_hat_t * r_t` (raw residual x raw regressor vector)

The papers (Wall & Amemiya 2003 Eq 13, 2000 Eq 28) define the exact version:
- `ell_t(alpha_hat) = m_hat_t - M_hat_t * alpha_hat`

where M_hat_t and m_hat_t are individual-level **bias-corrected** contributions (same recursion as point estimation, but per-observation instead of averaged).

### Files to modify: `R/est_2smm.R` and `R/modsem_2smm.R`

### Implementation Steps (in order):

1. **Add `make_error_moment_fn(W, eps_cumulants)` factory** in `est_2smm.R`
   - Creates memoized error_moment_fn closures from W and eps_cumulants
   - Replaces the inline closure in modsem_2smm.R (lines ~140-155)
   - Needed to create perturbed error_moment_fn during C_hat computation

2. **Add `compute_bc_individual(r, f_hat, error_moment_fn, memo)`** in `est_2smm.R`
   - Same recursion as existing `unbiased_moment()` (line 284) but returns an **n-vector** (one value per observation) instead of a scalar mean
   - `bc_t(r) = raw_product_t(r) - sum_{s>0, |s|>=2} C(r,s) * mu_e(s) * bc_t(r-s)`
   - For sum(r)=0: rep(1,n). For sum(r)=1: raw product (no correction since E[e]=0)
   - Memoized on multi-index key

3. **Add `compute_bc_ell(f_hat, alpha_hat, reg_multiidx, eta_name, factor_names, error_moment_fn)`** in `est_2smm.R`
   - Returns n x n_reg matrix
   - For each j: `ell_t[j] = bc_t(r_j + r_eta) - sum_k alpha_k * bc_t(r_j + r_k)`
   - This is the exact paper definition: ell_i = m_hat_i - M_hat_i * alpha

4. **Modify `compute_ell_contributions()`** (line 751 in est_2smm.R)
   - Add `error_moment_fn` parameter (default NULL)
   - If provided: call `compute_bc_ell` (bias-corrected ell)
   - If NULL: fall back to current raw ell (backward compatible)

5. **Modify `compute2SMMSE()`** (line 547 in est_2smm.R)
   - Add `ell` parameter (default NULL)
   - If provided: use it for Omega_star; otherwise compute raw ell internally
   - Raw residuals still computed for variance estimation output

6. **Modify `reconstruct_stage1_from_perturbed()`** (line 703 in est_2smm.R)
   - Also return `Theta` in the list (needed to update eps_cumulants[[2]] under perturbation)

7. **Modify `compute_C_hat()`** (line 787 in est_2smm.R)
   - Add `eps_cumulants` parameter
   - For each perturbation: build perturbed eps_cumulants (update [[2]] = diag(Theta_perturbed), keep [[3]], [[4]] fixed), create perturbed error_moment_fn via factory, pass to `compute_ell_contributions`
   - This captures BOTH paths from the chain rule (Eq 14): f_hat path + Sigma_ee path

8. **Wire in `modsem_2smm()`** (modsem_2smm.R)
   - Use `make_error_moment_fn` factory for the original error_moment_fn
   - After stage2, compute bias-corrected ell via `compute_bc_ell` and pass to `compute2SMMSE`
   - Pass `eps_cumulants` to `compute_C_hat`

### Key Insight
The individual-level bias-corrected contribution `bc_t(r)` uses the SAME recursion as `unbiased_moment(r)`, but per-observation:
- `unbiased_moment(r)` returns a scalar (mean across observations)
- `compute_bc_individual(r)` returns an n-vector (one per observation)

The error moment mu_e(s) terms are **constants** (population-level), so the only difference is that raw products are individual-level vectors instead of their means.

### Verification
1. `devtools::test()` — all existing tests pass
2. Run `simulations/quick_se_check.R` — SE/SD should remain ~1.0, coverage ~0.95
3. Results may improve slightly for cases where the approximation was noticeable
