### simulate_shared_controls.R
### Copyright 2022 Nicholas J. Seewald, PhD

#' Simulate Stacked DiD Cohorts with Shared Control Individuals
#'
#' Generate, summarize, and analyze a two-cohort stacked
#' difference-in-differences study.
#'
#' @param Tpre Number of measurement occasions in the pre-treatment period.
#'   Coerced to integer.
#' @param Tpost Number of measurement occasions in the post-treatment period.
#'   Coerced to integer.
#' @param Delta Number of measurement occasions between cohort study period
#'   start times. Coerced to integer.
#' @param nctrlstates Number of control states; one number common to both
#'   cohorts. Coerced to integer.
#' @param nperstate_cohort1 Number of individuals in each state in cohort 1.
#'   Recycled to length (`nctrlstates` + 1) such that each element specifies
#'   sample size for each state in the order (cohort 1's treated state, control
#'   state 1, control state 2, ...)
#' @param nperstate_cohort2 Number of individuals in each state in cohort 2; see
#'   `nperstate_cohort1` above.
#' @param nshared Number of shared control individuals per control state that
#'   contribute to analysis for both cohorts, of the full state sample size
#'   provided in the `nperstate` variables. Recycled to length `nctrlstates` and
#'   assumed to be given in the order (control state 1, control state 2, ...)
#' @param rho Exchangeable within-person correlation. If a single number,
#'   assumed constant across all states; else, a vector of length 2 +
#'   `nctrlstates` in the order (cohort 1 treated, cohort 2 treated, control
#'   states).
#' @param phi Within-period correlation, i.e., correlation between two
#'   observations from different people in the same state at the same time.
#'   Recycled to length 2 + `nctrlstates` in the order (cohort 1 treated, cohort
#'   2 treated, control states).
#' @param psi Within-state correlation, i.e., correlation between two
#'   observations from different people in the same state at different times.
#'   Recycled to length 2 + `nctrlstates` in the order (cohort 1 treated, cohort
#'   2 treated, control states).
#' @param tx_intshift_cohort1 Numeric. Intercept shift at time of treatment
#'   \eqn{t^*} for the treated state in cohort 1.
#' @param tx_intshift_cohort2 Numeric. Intercept shift at time of treatment
#'   \eqn{t^*} for the treated state in cohort 2.
#' @param state_fe_mean Mean of state fixed effects. See Details.
#' @param error_sd Standard deviation of additive error in data generative
#'   model, assumed common across person-state-time. See Details.
#' @param secular_trend Function of one argument (time) defining a secular time
#'   trend common to all states in the analysis. Defaults to 0.
#' @param ARdecay Autoregressive parameter for the correlation structure of the
#'   unit-time random effects in the data-generative model. Must be between 0
#'   and 1; 1 indicates no autoregressive decay, 0 implies unit-time random
#'   effects are pairwise independent.
#'
#' @details Generates, summarizes, and analyzes a `data.frame` for a two-cohort
#'   stacked difference-in-differences analysis with shared control individuals.
#'   The `data.frame` contains individual-level data from two treated states and
#'   `nctrlstates` control states that are common to both treated states'
#'   analyses.
#'
#'   The generative model for the outcome observed from individual \eqn{i} at
#'   time \eqn{} unit/state \eqn{u} is \deqn{Y_{uit} = \beta_0 + \beta_{1,t} +
#'   \beta_2 A_{u t} + b_{i} + b_{u} + b_{u t} + \epsilon_{uit},} where
#'   \eqn{\beta_0} corresponds to `state_fe_mean`, \eqn{\beta_{1,t} =
#'   \beta_1(t)} is `secular_trend(`\eqn{t}`)`; \eqn{\beta_2} is
#'   `tx_intshift_cohortK`; \eqn{A_{ut}} is an indicator for whether unit
#'   \eqn{u} is treated at time \eqn{t};  \eqn{b_i}, \eqn{b_u}, and \eqn{b_{ut}}
#'   are random effects for individual, unit, and unit-time, respectively; and
#'   \eqn{\epsilon_{uit}} is residual error with standard deviation `error_sd`.
#'   In the above, `K` is either 1 or 2 depending on the cohort for which data
#'   is being generated. Note that state fixed effects are *not*
#'   user-specified: they are defined as \eqn{\beta_0 + b_u}, with \eqn{b_u \sim
#'   N(0, \sigma^2_{bu})}. Also, note that the treatment effect is *constant*
#'   over time.
#'   
#'   The analytic 
#'
#' @importFrom checkmate makeAssertCollection assert_integerish assert_numeric
#'   assert_function reportAssertions assert_number assert_count
#'
#' @return A list
#' @export
#'
#' @examples
simulate_shared_controls <- function(Tpre, Tpost, Delta, 
                                     nctrlstates, 
                                     nperstate_cohort1,
                                     nperstate_cohort2 = nperstate_cohort1,
                                     nshared,
                                     rho = rep(0, 2 + nctrlstates),
                                     phi = rep(0, 2 + nctrlstates),
                                     psi = rep(0, 2 + nctrlstates),
                                     tx_intshift_cohort1,
                                     tx_intshift_cohort2 = tx_intshift_cohort1,
                                     secular_trend = function(x) 0,
                                     state_fe_mean = 0,
                                     error_sd = 1,
                                     ARdecay = 1) {
  # CHECK FOR WELL-FORMED INPUT
  checks <- makeAssertCollection()
  
  # Check design arguments: must be single integers >= 1
  assert_count(Tpre, positive = T, add = checks)
  assert_count(Tpost, positive = T, add = checks)
  assert_count(Delta, positive = T, add = checks)
  
  # Check sample size arguments: must be integers >= 1
  assert_integerish(nctrlstates, lower = 1, add = checks)
  assert_integerish(nperstate_cohort1, lower = 1, max.len = nctrlstates + 1,
                    add = checks)
  assert_integerish(nperstate_cohort2, lower = 1, max.len = nctrlstates + 1,
                    add = checks)
  assert_integerish(nshared, lower = 1, max.len = nctrlstates + 1, add = checks)
  
  # Recycle sample size vectors if needed
  if ((length(nperstate_cohort1) %% (nctrlstates + 1)) > 1)
    warning("nperstate_cohort1 is not length 1 or (nctrlstates + 1). Recycling.")
  nperstate_cohort1 <- rep_len(nperstate_cohort1, nctrlstates + 1)
  
  if ((length(nperstate_cohort2) %% (nctrlstates + 1)) > 1)
    warning("nperstate_cohort2 is not length 1 or (nctrlstates + 1). Recycling.")
  nperstate_cohort2 <- rep_len(nperstate_cohort2, nctrlstates + 1)
  
  if ((length(nshared) %% nctrlstates) > 1)
    warning("nshared is not length 1 or nctrlstates. Recycling.")
  nshared <- rep_len(nshared, nctrlstates)
  
  # Make sure nshared isn't larger than number of individuals in that state
  assert_numeric(pmin(nperstate_cohort1[-1], nperstate_cohort2[-1]) - nshared,
                 lower = 0, add = checks)
  
  # Check correlation arguments
  assert_numeric(rho, upper = 1, lower = -1, max.len = nctrlstates + 2, 
                 any.missing = F, all.missing = F, add = checks)
  assert_numeric(phi, upper = 1, lower = -1, max.len = nctrlstates + 2, 
                 any.missing = F, all.missing = F, add = checks)
  assert_numeric(psi, upper = 1, lower = -1, max.len = nctrlstates + 2, 
                 any.missing = F, all.missing = F, add = checks)
  
  # Recycle correlation vectors if needed
  if ((length(rho) %% (nctrlstates + 2)) > 1)
    warning("rho is not length 1 or (nctrlstates + 2). Recycling.")
  rho <- rep_len(rho, nctrlstates + 2)
  
  if ((length(phi) %% (nctrlstates + 2)) > 1)
    warning("phi is not length 1 or (nctrlstates + 2). Recycling.")
  phi <- rep_len(phi, nctrlstates + 2)
  
  if ((length(psi) %% (nctrlstates + 2)) > 1)
    warning("psi is not length 1 or (nctrlstates + 2). Recycling.")
  psi <- rep_len(psi, nctrlstates + 2)
  
  if ((length(nshared) %% nctrlstates) > 1)
    warning("nshared is not length 1 or nctrlstates. Recycling.")
  nshared <- rep_len(nshared, nctrlstates)
  
  # Require psi <= phi <= rho
  assert_numeric(rho - phi, lower = 0, add = checks)
  assert_numeric(phi - psi, lower = 0, add = checks)
  
  # Check remaining arguments
  assert_numeric(ARdecay, lower = 0, upper = 1, any.missing = F, add = checks)
  assert_function(secular_trend, nargs = 1, add = checks)
  assert_number(state_fe_mean, finite = T, add = checks)
  assert_number(error_sd, lower = 0, add = checks)
  error_sd <- rep_len(error_sd, length.out = nctrlstates + 2)
  
  # Report collected input errors
  reportAssertions(checks)
  
  # Convert to integer
  # Tpre <- as.integer(Tpre)
  # Tpost <- as.integer(Tpost)
  # Delta <- as.integer(Delta)
  # nperstate_cohort1 <- as.integer(nperstate_cohort1)
  # nperstate_cohort2 <- as.integer(nperstate_cohort2)
  # nshared <- as.integer(nshared)
  
  # Adjust nperstate in cohort 2 to generate different IDs than in cohort 1 (IDs
  # are a hash of a person index and state name, so given a control state
  # generate_cohort will by default return the same control IDs every time. Here
  # we generate more IDs than necessary so we can delete the reused ones that
  # aren't shared controls)
  # nperstate_cohort2[-1] <- 2 * nperstate_cohort2[-1] - nshared
  
  # Generate early cohort
  c1 <- generate_cohort(nperstate = nperstate_cohort1,
                        nctrlstates = nctrlstates,
                        Tobs = Tpre + Tpost,
                        Tpost = Tpost, 
                        txstatename = "tx01",
                        cohort_name = "cohort1",
                        start_time = 0)
  
  # Select IDs in each state in cohort 1 to also appear in cohort 2
  shared_ids <- c(sapply(1:nctrlstates, \(i) {
    sample(
      unique(
        c1$cohort$id[
          c1$cohort$state == paste0("ctrl", 
                                    formatC(i, 
                                            width = nchar(nctrlstates),
                                            flag = "0"))]),
      size = nshared[i], replace = F)
  }))
  
  # Generate late cohort: build skeleton starting at time Delta
  c2 <- generate_cohort(nperstate = nperstate_cohort2,
                        nctrlstates = nctrlstates,
                        Tobs = Tpre + Tpost,
                        Tpost = Tpost, 
                        txstatename = "tx02",
                        cohort_name = "cohort2",
                        start_time = Delta)
  
  # Generate late cohort: sample IDs to remove from control states and replace
  # with IDs from early cohort
  removed_ids <- c(sapply(1:nctrlstates, \(i) {
    sample(
      unique(
        c2$cohort$id[
          c2$cohort$state == paste0("ctrl", 
                                    formatC(i, 
                                            width = nchar(nctrlstates),
                                            flag = "0"))
        ]),
      size = nshared[i], replace = F)
  }))
  
  # Generate late cohort: replace removed IDs with shared IDs
  c2$cohort$id[c2$cohort$id %in% removed_ids] <- 
    rep(shared_ids, each = Tpre + Tpost)

  # Vector of state names
  states <- c("tx01", "tx02", 
              paste0(c(rep("ctrl", nctrlstates)),
                     formatC(1:nctrlstates, width = nchar(nctrlstates),
                             flag = "0")))

  # c1$cohort <- c1$cohort[order(c1$cohort$state, c1$cohort$id, c1$cohort$time), ]
  # c2$cohort <- c2$cohort[order(c2$cohort$state, c2$cohort$id, c2$cohort$time), ]

  c1$cohort <- transform(c1$cohort, 
                         "Ymean" = state_fe_mean + secular_trend(time) + 
                           tx_intshift_cohort1 * treated)
  c2$cohort <- transform(c2$cohort, 
                         "Ymean" = state_fe_mean + secular_trend(time) +
                           tx_intshift_cohort2 * treated)

  c1$cohort$Y <- c2$cohort$Y <- NA

  #### RANDOM EFFECTS TO INDUCE CORRELATION (Kasza et al. 2019)
  
  # Get random effect variances
  state_varcors <- compute_raneff_vars(rho, phi, psi, error_sd, c1, c2)

  # Get vector of all times in either cohort
  fullTimes <- seq(min(c1$cohort$time), max(c2$cohort$time))
  
  # Parametrize R matrices (correlation btwn state-time random effects): this is
  # a T-by-T matrix (T = Tpre + Tpost) where the (t,s) element is the
  # correlation between the random effect at time t and the random effect of
  # time s. It's a multiplier on the within-period ICC.
  R <- structure(
    lapply(state_varcors$r, \(x) {
      # For each state, make a matrix of r = psi/phi with 1's on diagonal
      temp <- matrix(x, nrow = length(fullTimes), ncol = length(fullTimes))
      diag(temp) <- 1
      # If autoregressive correlation specified, make an AR1(ARdecay) matrix
      if (ARdecay != 1) {
        decay <- matrix(1, nrow = length(fullTimes), ncol = length(fullTimes))
        for (i in 1:length(fullTimes)) {
          for (j in 1:i) {
            decay[i, j] <- decay[j, i] <- ARdecay^(abs(i - j))
          }
        }
        # multiply the decay matrix by the temp r matrix from above to let r
        # decay over time
        return(temp * decay)
      } else return(temp)
    }),
    names = states)
  
  # Generate "cluster-period" (i.e., state-time) random effects
  CP <- do.call(
    rbind,
    lapply(1:length(R), \(i) {
      CPu <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(R[[i]])),
                           Sigma = state_varcors$s2_sttime[i] * R[[i]])
      data.frame("state" = states[i], "time" = fullTimes,
                 "raneff_statetime" = CPu)
    }))

  # Make data.frame of random effects for state & person
  raneff <- merge(c1$cohort[, c("id", "state")], c2$cohort[, c("id", "state")],
                  all = T)
  raneff <- raneff[!duplicated(raneff), ]
  
  
  raneff <- do.call(rbind, lapply(split.data.frame(raneff, raneff$state), \(x) {
    x$raneff_id <- rnorm(nrow(x), mean = 0,
                      sd = sqrt(state_varcors$s2_id[state_varcors$state ==
                                                      unique(x$state)]))
    x
  }))
  
  outcome_error <-
    expand.grid(id = union(c1$cohort$id, c2$cohort$id),
                time = union(c1$cohort$time, c2$cohort$time))

  outcome_error$error <- rnorm(n = nrow(outcome_error),
                               mean = 0, sd = error_sd)

  # merge state-time random effects into cohort data.frames
  c1$cohort <- merge(c1$cohort, CP, by = c("state", "time"), all.x = T) |>
    merge(raneff, by = c("state", "id"), all.x = T)
  c2$cohort <- merge(c2$cohort, CP, by = c("state", "time"), all.x = T)|>
    merge(raneff, by = c("state", "id"), all.x = T)
  
  c1$cohort <- merge(c1$cohort, outcome_error, all.x = T, by = c("id", "time"))
  c2$cohort <- merge(c2$cohort, outcome_error, all.x = T, by = c("id", "time"))

  # Add random effects and noise to Ymean
  c1$cohort$Y <- with(c1$cohort,
                      Ymean + raneff_statetime + raneff_id + error)
  c2$cohort$Y <- with(c2$cohort,
                      Ymean + raneff_statetime + raneff_id + error)
  
  c1c2t <- time_windows(unique(c1$cohort$time), unique(c2$cohort$time),
                        c1tStar = c1$tStar, c2tStar = c2$tStar)
  periods <- do.call(rbind, c1c2t)
  periods <- data.frame("time" = periods[, 1], "period" = rownames(periods),
                        "shared" = "shared")

  c1$cohort$shared <- ifelse(c1$cohort$id %in% shared_ids, "shared", "disjoint")
  c2$cohort$shared <- ifelse(c2$cohort$id %in% shared_ids, "shared", "disjoint")

  c1$cohort <- merge(c1$cohort, periods, all.x = T, by = c("time", "shared"))
  c2$cohort <- merge(c2$cohort, periods, all.x = T, by = c("time", "shared"))

  c1$cohort <- c1$cohort |>
    transform("period" = ifelse(shared == "shared",
                                period, ifelse(post, "post", "pre")))
  c2$cohort <- c2$cohort |>
    transform("period" = ifelse(shared == "shared",
                                period, ifelse(post, "post", "pre")))

  # Expand cohort data frames into lists of data frames, one per subgroup of
  # interest (i.e., shared control individuals in the post/post window, etc.)
  c1_expand <- split.data.frame(c1$cohort, f = ~ cohort + state + shared + period)
  c1_expand <- c1_expand[sapply(c1_expand, \(y) nrow(y) != 0)]
  c2_expand <- split.data.frame(c2$cohort, f = ~ cohort + state + shared + period)
  c2_expand <- c2_expand[sapply(c2_expand, \(y) nrow(y) != 0)]
  
  c1disj   <- subset(c1$cohort, !(id %in% shared_ids) & state != "tx01")
  c1shared <- subset(c1$cohort, id %in% shared_ids & state != "tx01")
  c2disj   <- subset(c2$cohort, !(id %in% shared_ids) & state != "tx02")
  c2shared <- subset(c2$cohort, id %in% shared_ids & state != "tx02")
  
  # Linear Models
  c1_lm <- lm(Y ~ -1 + state + as.factor(time) + treated, data = c1$cohort)
  c2_lm <- lm(Y ~ -1 + state + as.factor(time) + treated, data = c2$cohort)
  
  # True overall treatment effect
  truth <- mean(c(tx_intshift_cohort1, tx_intshift_cohort2))
  
  component_sums <- sapply(c(c1_expand, c2_expand), \(x) sum(x$Y, na.rm = T))
  
  component_sums_old <- data.frame(
    "c1_tx_pre"    = with(c1$cohort, sum(Y[state == "tx01" & !post])),
    "c1_tx_post"   = with(c1$cohort, sum(Y[state == "tx01" & post])),
    "c1_disj_pre"  = with(c1disj,    sum(Y[!post])),
    "c1_disj_post" = with(c1disj,    sum(Y[post])),
    "c1_pre_dot"   = with(c1shared,  sum(Y[time %in% c1c2t$pre_dot])),
    "c1_pre_pre"   = with(c1shared,  sum(Y[time %in% c1c2t$pre_pre])),
    "c1_post_pre"  = with(c1shared,  sum(Y[time %in% c1c2t$post_pre])),
    "c1_post_post" = with(c1shared,  sum(Y[time %in% c1c2t$post_post])),
    "c1_post_dot"  = with(c1shared,  sum(Y[time %in% c1c2t$post_dot])),
    "c2_tx_pre"    = with(c2$cohort, sum(Y[state == "tx02" & !post])),
    "c2_tx_post"   = with(c2$cohort, sum(Y[state == "tx02" & post])),
    "c2_disj_pre"  = with(c2disj,    sum(Y[!post])),
    "c2_disj_post" = with(c2disj,    sum(Y[post])),
    "c2_pre_pre"   = with(c2shared,  sum(Y[time %in% c1c2t$pre_pre])),
    "c2_post_pre"  = with(c2shared,  sum(Y[time %in% c1c2t$post_pre])),
    "c2_post_post" = with(c2shared,  sum(Y[time %in% c1c2t$post_post])),
    "c2_dot_pre"   = with(c2shared,  sum(Y[time %in% c1c2t$dot_pre])),
    "c2_dot_post"  = with(c2shared,  sum(Y[time %in% c1c2t$dot_post]))
  )
  
  component_means <- data.frame(
    "c1_tx_pre"    = with(c1$cohort, mean(Y[state == "tx01" & !post])),
    "c1_tx_post"   = with(c1$cohort, mean(Y[state == "tx01" & post])),
    "c1_ctrl_pre"  = with(c1$cohort, mean(Y[state != "tx01" & !post])),
    "c1_ctrl_post" = with(c1$cohort, mean(Y[state != "tx01" & post])),
    "c2_tx_pre"    = with(c2$cohort, mean(Y[state == "tx02" & !post])),
    "c2_tx_post"   = with(c2$cohort, mean(Y[state == "tx02" & post])),
    "c2_ctrl_pre"  = with(c2$cohort, mean(Y[state != "tx02" & !post])),
    "c2_ctrl_post" = with(c2$cohort, mean(Y[state != "tx02" & post]))
  )
  
  did_estimates <- data.frame(
    "c1Est_lm" = coef(c1_lm)["treatedTRUE"],
    "c2Est_lm" = coef(c2_lm)["treatedTRUE"],
    "c1Est_emp" = with(component_means, (c1_tx_post - c1_tx_pre) - 
                      (c1_ctrl_post - c1_ctrl_pre)),
    "c2Est_emp" = with(component_means, (c2_tx_post - c2_tx_pre) - 
                      (c2_ctrl_post - c2_ctrl_pre))
  )
  
  did_ses <- data.frame(
    "c1SE_lm" = summary(c1_lm)$coefficients["treatedTRUE", 2],
    "c2SE_lm" = summary(c2_lm)$coefficients["treatedTRUE", 2],
    "c1SE_idAdj" = 
      sqrt(sandwich::vcovCL(c1_lm,
                            cluster = ~ id)["treatedTRUE", "treatedTRUE"]),
    "c2SE_idAdj" = 
      sqrt(sandwich::vcovCL(c2_lm,
                            cluster = ~ id)["treatedTRUE", "treatedTRUE"]),
    "c1SE_stAdj" = 
      sqrt(sandwich::vcovCL(c1_lm,
                            cluster = ~ state)["treatedTRUE", "treatedTRUE"]),
    "c2SE_stAdj" = 
      sqrt(sandwich::vcovCL(c2_lm,
                            cluster = ~ state)["treatedTRUE", "treatedTRUE"]),
    "c1SE_analytic" = sqrt(
      varDD(ntx = nperstate_cohort1[1], 
            nctrl = nperstate_cohort1[-1],
            Tpre = Tpre, Tpost = Tpost,
            rho = with(state_varcors, rho[state != "tx02"]), 
            phi = with(state_varcors, phi[state != "tx02"]),
            psi = with(state_varcors, psi[state != "tx02"]),
            sigma_s = with(state_varcors, 
                           sqrt(outcome_var[state != "tx02"])))),
    "c2SE_analytic" = sqrt(
      varDD(ntx = nperstate_cohort2[1], 
            nctrl = nperstate_cohort2[-1],
            Tpre = Tpre, Tpost = Tpost,
            rho = with(state_varcors, rho[state != "tx01"]), 
            phi = with(state_varcors, phi[state != "tx01"]),
            psi = with(state_varcors, psi[state != "tx01"]),
            sigma_s = with(state_varcors, 
                           sqrt(outcome_var[state != "tx01"])))),
    "c1c2cor_analytic" = 
      corDD(Ntx_cohort1 = nperstate_cohort1[1],
            Ntx_cohort2 = nperstate_cohort2[1],
            Ndisj_cohort1 = nperstate_cohort1[-1] - nshared,
            Ndisj_cohort2 = nperstate_cohort2[-1] - nshared,
            Nshared = nshared,
            Tpre = Tpre, Tpost = Tpost, Delta = Delta,
            rho = rho, phi = phi, psi = psi,
            outcomeSD = sqrt(state_varcors$outcome_var))
  )
  
  V_analytic <- diag(c(did_ses$c1SE_analytic, did_ses$c2SE_analytic)) %*% 
    cormat(rho = did_ses$c1c2cor_analytic, p = 2, corstr = "exch") %*%
    diag(c(did_ses$c1SE_analytic, did_ses$c2SE_analytic))
  
  # Aggregate via inverse-variance weighting using unadjusted SEs from lm()
  ivw_lmSE <- data.frame(
    "ivw_lmSE_est" = 
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)),
        w = with(did_ses, c(1/c1SE_lm^2, 1/c2SE_lm^2))),
    "ivw_lmSE_uncorSE" = 
      1 / sqrt(with(did_ses, c1SE_lm^(-2) + c2SE_lm^(-2))),
    "ivw_lmSE_corSE" =
      with(did_ses, c(c1SE_lm^(-2), c2SE_lm^(-2)) %*% V_analytic %*%
             c(c1SE_lm^(-2), c2SE_lm^(-2)))^(-1/2)
  ) |>
    transform(
      "ivw_lmSE_uncorLB" = ivw_lmSE_est - qnorm(.975) * ivw_lmSE_uncorSE,
      "ivw_lmSE_uncorUB" = ivw_lmSE_est + qnorm(.975) * ivw_lmSE_uncorSE,
      "ivw_lmSE_corLB"   = ivw_lmSE_est - qnorm(.975) * ivw_lmSE_corSE,
      "ivw_lmSE_corLB"   = ivw_lmSE_est - qnorm(.975) * ivw_lmSE_corSE,
    ) |> 
    transform(
      "ivw_lmSE_uncorCovg" = truth >= ivw_lmSE_uncorLB & truth <= ivw_lmSE_uncorUB,
      "ivw_lmSE_corCovg".  = truth >= ivw_lmSE_corLB.  & truth <= ivw_lmSE_corUB
    )
  
  # Aggregate via inverse-variance weighting using cluster-adjusted SEs with
  # clustering at the individual level
  ivw_idAdjSE <- data.frame(
    "ivw_idAdjSE_est" = 
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)),
        w = with(did_ses, c(1/c1SE_idAdj^2, 1/c2SE_idAdj^2))),
    "ivw_idAdjSE_uncorSE" = 
      1 / sqrt(with(did_ses, c1SE_idAdj^(-2) + c2SE_idAdj^(-2))),
    "ivw_idAdjSE_corSE" =
      with(did_ses, c(c1SE_idAdj^(-2), c2SE_idAdj^(-2)) %*% V_analytic %*%
             c(c1SE_idAdj^(-2), c2SE_idAdj^(-2)))^(-1/2)
  ) |>
    transform(
      "ivw_lmSE_uncorLB" = ivw_lmSE_est - qnorm(.975) * ivw_lmSE_uncorSE,
      "ivw_lmSE_uncorUB" = ivw_lmSE_est + qnorm(.975) * ivw_lmSE_uncorSE,
      "ivw_lmSE_corLB"   = ivw_lmSE_est - qnorm(.975) * ivw_lmSE_corSE,
      "ivw_lmSE_corLB"   = ivw_lmSE_est - qnorm(.975) * ivw_lmSE_corSE,
    ) |> 
    transform(
      "ivw_lmSE_uncorCovg" = truth >= ivw_lmSE_uncorLB & truth <= ivw_lmSE_uncorUB,
      "ivw_lmSE_corCovg"  = truth >= ivw_lmSE_corLB.  & truth <= ivw_lmSE_corUB
    )
  
  ivw_estimates <- data.frame(
    "truth" = truth,
    "ivw_lmSE_est" = 
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)),
        w = with(did_ses, c(1/c1SE_lm^2, 1/c2SE_lm^2))),
    "ivw_lmSE_uncorSE" = 
      1 / sqrt(with(did_ses, c1SE_lm^(-2) + c2SE_lm^(-2))),
    "ivw_idAdjSE_est" =
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)),
        w = with(did_ses, c(1/c1SE_idAdj^2, 1/c2SE_idAdj^2))),
    "ivw_stAdjSE_est" = 
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)),
        w = with(did_ses, c(1/c1SE_stAdj^2, 1/c2SE_stAdj^2))),
    "ivw_analyticSE_est" = 
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)),
        w = with(did_ses, c(1/c1SE_analytic^2, 1/c2SE_analytic^2))),
    
    "ivw_idAdjSE_uncorSE" = 1 / sqrt(with(did_ses, c1SE_idAdj^(-2) + c2SE_idAdj^(-2))),
    "ivw_stAdjSE_uncorSE" = 1 / sqrt(with(did_ses, c1SE_stAdj^(-2) + c2SE_stAdj^(-2))),
    "ivw_analyticSE_uncorSE" = 1 / sqrt(with(did_ses, 
                                  c1SE_analytic^(-2) + c2SE_analytic^(-2))),
    "ivw_idAdj_SE"
    ) |>
    transform(
              "ivw_idAdjSE_lb"    = ivw_idAdjSE_est - qnorm(.975) * ivw_idAdjSE_uncorSE,
              "ivw_idAdjSE_ub"    = ivw_idAdjSE_est + qnorm(.975) * ivw_idAdjSE_uncorSE,
              "ivw_stAdjSE_lb"    = ivw_stAdjSE_est - qnorm(.975) * ivw_stAdjSE_uncorSE,
              "ivw_stAdjSE_ub"    = ivw_stAdjSE_est + qnorm(.975) * ivw_stAdjSE_uncorSE,
              "ivw_analyticSE_lb" = ivw_analyticSE_est - qnorm(.975) * ivw_analyticSE_uncorSE,
              "ivw_analyticSE_ub" = ivw_analyticSE_est + qnorm(.975) * ivw_analyticSE_uncorSE
    ) |> 
    transform("ivw_lmSE_covg"       = truth >= ivw_lmSE_lb & truth <= ivw_lmSE_ub,
              "ivw_idAdjSE_covg"    = truth >= ivw_idAdjSE_lb & truth <= ivw_idAdjSE_ub,
              "ivw_stAdtSE_covg"    = truth >= ivw_stAdjSE_lb & truth <= ivw_stAdjSE_ub,
              "ivw_analyticSE_covg" = truth >= ivw_analyticSE_lb & truth <= ivw_analyticSE_ub)
  
  gls_estimates <- data.frame(
    "gls_est" =
      weighted.mean(
        x = with(did_estimates, c(c1Est_lm, c2Est_lm)), 
        w = c(V_analytic[2, 2] - V_analytic[2, 1],
              V_analytic[1, 1] - V_analytic[1, 2])
      ),
    "gls_SE" = sum(solve(V_analytic))^(-1/2)
  ) |> 
    transform(
      "gls_lb" = gls_est - qnorm(.975) * gls_SE,
      "gls_ub" = gls_est + qnorm(.975) * gls_SE
    ) |> 
    transform(
      "gls_covg" = truth >= gls_lb & truth <= gls_ub
    )
  
  out <- list("component_sums"  = component_sums,
              "component_sums_old" = component_sums_old,
              "component_means" = component_means,
              "did_estimates"   = did_estimates,
              "did_ses" = did_ses,
              "ivw_estimates" = ivw_estimates,
              "gls_estimates" = gls_estimates)
  
  class(out) <- c("sharedCtrlsSim", class(out))
  
  return(out)
}

print.sharedCtrlsSim <- function(x) {
  cat("\t Simulated DiD Results with Shared Control Individuals\n\n")
  cat()
}
