### simulate_shared_controls.R
### Copyright 2022 Nicholas J. Seewald, PhD

#' Simulate Stacked DiD Cohorts with Shared Control Individuals
#'
#' @param Tpre Number of measurement occasions in the pre-treatment period
#' @param Tpost Number of measurement occasions in the post-treatment period
#' @param Delta Number of measurement occasions between cohort study period
#'   start times
#' @param nctrlstates Number of control states in each cohort
#' @param nperstate_cohort1 Number of individuals in each state in cohort 1.
#'   Optionally, a vector of length 2 + nctrlstates that specifies different
#'   numbers for each state in the order (cohort 1 treated, cohort 2 treated,
#'   control states). Coerced to integer.
#' @param nperstate_cohort2 Number of individuals in each state in cohort 2; see
#'   above.
#' @param nshared Number of shared control individuals
#' @param rho Exchangeable within-person correlation. If a single number,
#'   assumed constant across all states; else, a vector of length 2 +
#'   nctrlstates in the order (cohort 1 treated, cohort 2 treated, control
#'   states).
#' @param phi Within-period correlation, i.e., correlation between two
#'   observations from different people in the same state at the same time. If a
#'   single number, assumed constant across all states; else, a vector of length
#'   2 + nctrlstates in the order (cohort 1 treated, cohort 2 treated, control
#'   states).
#' @param psi Within-state correlation, i.e., correlation between two
#'   observations from different people in the same state at different times. If
#'   a single number, assumed constant across all states; else, a vector of
#'   length 2 + nctrlstates in the order (cohort 1 treated, cohort 2 treated,
#'   control states).
#' @param tx_intshift_cohort1 Numeric. Intercept shift at time of treatment
#'   \eqn{t^*} for the treated state in cohort 1.
#' @param tx_intshift_cohort2 Numeric. Intercept shift at time of treatment
#'   \eqn{t^*} for the treated state in cohort 2.
#' @param secular_cohort1 Function of one argument defining a secular time trend
#'   common to all states in cohort 1. Defaults to 0.
#' @param secular_cohort2 Function of one argument defining a secular time trend
#'   common to all states in cohort 2. Defaults to 0.
#' @param state_fe_mean Mean of state fixed effects. See Details.
#' @param error_sd Standard deviation of additive error in data generative
#'   model. See Details.
#' @param verbose Logical. Should additional messages be printed?
#'
#' @details Generates and analyzes a `data.frame` for stacked
#'   difference-in-differences analysis with shared control individuals. The
#'   `data.frame` contains individual-level data from two treated states and
#'   `nctrlstates` control states that are common to both treated states' analyses. 
#'
#' @return
#' @export
#'
#' @examples
simulate_shared_controls <- function(Tpre, Tpost, Delta, 
                                     nctrlstates, 
                                     nperstate_cohort1,
                                     nperstate_cohort2 = nperstate_cohort1,
                                     nshared,
                                     rho = rep(0, 2 + nctrlstates),
                                     phi = rep(0, 2 + nctrlstates), ## CHECK LENGTH
                                     psi = rep(0, 2 + nctrlstates), ## CHECK LENGTH
                                     tx_intshift_cohort1,
                                     tx_intshift_cohort2 = tx_intshift_cohort1,
                                     secular_trend = function(x) 0,
                                     state_fe_mean = 0,
                                     error_sd = rep(1, 2 + nctrlstates),
                                     verbose = FALSE,
                                     ARdecay = 1) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(rho, add = checks)
  checkmate::assert_numeric(phi, add = checks)
  checkmate::assert_numeric(psi, add = checks)
  checkmate::reportAssertions(checks)
  
  rho <- expand_vec(rho, 2 + nctrlstates)
  phi <- expand_vec(phi, 2 + nctrlstates)
  psi <- expand_vec(psi, 2 + nctrlstates)
  error_sd <- expand_vec(error_sd, 2 + nctrlstates)

  
  if (length(nperstate_cohort1) == 1) {
    nperstate_cohort1 <- rep_len(nperstate_cohort1, nctrlstates + 1)
  } else if (length(nperstate_cohort1) != nctrlstates + 1)
    stop("nperstate_cohort1 must have length 1 or (nctrlstates + 1)")
  
  if (length(nperstate_cohort2) == 1) {
    nperstate_cohort2 <- rep_len(nperstate_cohort2, nctrlstates + 1)
  } else if (length(nperstate_cohort2) != nctrlstates + 1)
    stop("nperstate_cohort2 must have length 1 or (nctrlstates + 1)")
  
  if (length(nshared) == 1) {
    nshared <- rep(nshared, nctrlstates)
  } else if (length(nshared) != nctrlstates) {
    stop("nshared must have length 1 or the same length as nctrlstates")
  }
  if (max(nshared) > min(c(nperstate_cohort1, nperstate_cohort2)))
    stop ("nshared is too big relative to number of people per state")
  if (length(unique(error_sd)) != 1)
    stop("Different error variances across states is currently unsupported.")
  
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
  
  # Select IDs in each state to share across cohorts
  shared_ids <- c(sapply(1:nctrlstates, \(i) {
    sample(unique(
      c1$cohort$id[c1$cohort$state == paste0("ctrl", 
                                             formatC(i, 
                                                     width = nchar(nctrlstates),
                                                     flag = "0"))]),
      size = nshared[i], replace = F)
  }))
  
  # Generate late cohort: build skeleton
  c2 <- generate_cohort(nperstate = nperstate_cohort2,
                        nctrlstates = nctrlstates,
                        Tobs = Tpre + Tpost,
                        Tpost = Tpost, 
                        txstatename = "tx02",
                        cohort_name = "cohort2",
                        start_time = Delta)
  
  # Generate late cohort: sample IDs to remove from control states and replace
  # with IDs from early cohort
  killed_ids <- c(sapply(1:nctrlstates, \(i) {
    sample(unique(
      c2$cohort$id[c2$cohort$state == paste0("ctrl", 
                                             formatC(i, 
                                                     width = nchar(nctrlstates),
                                                     flag = "0"))]),
      size = nshared[i], replace = F)
  }))
  
  # Generate late cohort: replace killed IDs with shared IDs
  c2$cohort$id[c2$cohort$id %in% killed_ids] <- rep(shared_ids, each = Tpre + Tpost)
  
  ### TEST: make big covariance matrix and generate noise from there ###

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
  
  # Start by parametrizing R matrix (correlation btwn state-time random
  # effects): this is a T-by-T matrix (T = Tpre + Tpost) where the (t,s) element
  # is the correlation between the random effect at time t and the random effect
  # of time s. It's a multiplier on the within-period ICC. 

  state_varcors <- compute_raneff_vars(rho, phi, psi, error_sd, c1, c2)

  fullTimes <- seq(min(c1$cohort$time), max(c2$cohort$time))

  R <- structure(
    lapply(state_varcors$r, \(x) {
      temp <- matrix(x, nrow = length(fullTimes), ncol = length(fullTimes))
      diag(temp) <- 1
      if (ARdecay != 1) {
        decay <- matrix(1, nrow = length(fullTimes), ncol = length(fullTimes))
        for (i in 1:length(fullTimes)) {
          for (j in 1:i) {
            decay[i, j] <- decay[j, i] <- ARdecay^(abs(i - j))
          }
        }
        return(temp * decay)
      } else return(temp)
    }),
    names = states)

  CP <- do.call(rbind,
                lapply(1:length(R), \(i) {
                  CPu <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(R[[i]])),
                                       Sigma = state_varcors$s2_sttime[i] * R[[i]])
                  data.frame("state" = states[i], "time" = fullTimes,
                             "raneff_statetime" = CPu)
                }))

  raneff <- merge(c1$cohort[, c("id", "state")], c2$cohort[, c("id", "state")], all = T)
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

  # c1$cohort$ftime <- factor(c1$cohort$time)
  # c1$cohort$stime <- paste0(c1$cohort$state, c1$cohort$time)
  #
  # mod1 <- lme4::lmer(Y ~ ftime + treated + (1 | stime), data = c1$cohort)
  # mod2 <- lme4::lmer(Y ~ ftime + treated + (1 | stime) + (1 | id), data = c1$cohort)
  # mod3 <- lme4::lmer(Y ~ treated + (1 | ftime) + (1 | stime) + (1 | id), data = c1$cohort)

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
    "c1_lm" = coef(c1_lm)["treatedTRUE"],
    "c2_lm" = coef(c2_lm)["treatedTRUE"],
    "c1_emp" = with(component_means, (c1_tx_post - c1_tx_pre) - 
                      (c1_ctrl_post - c1_ctrl_pre)),
    "c2_emp" = with(component_means, (c2_tx_post - c2_tx_pre) - 
                      (c2_ctrl_post - c2_ctrl_pre))
  )
  
  did_ses <- data.frame(
    "c1_lm" = summary(c1_lm)$coefficients["treatedTRUE", 2],
    "c2_lm" = summary(c2_lm)$coefficients["treatedTRUE", 2],
    "c1_lm_idadj" = 
      sqrt(sandwich::vcovCL(c1_lm,
                            cluster = ~ id)["treatedTRUE", "treatedTRUE"]),
    "c2_lm_idadj" = 
      sqrt(sandwich::vcovCL(c2_lm,
                            cluster = ~ id)["treatedTRUE", "treatedTRUE"]),
    "c1_lm_stadj" = 
      sqrt(sandwich::vcovCL(c1_lm,
                            cluster = ~ state)["treatedTRUE", "treatedTRUE"]),
    "c2_lm_stadj" = 
      sqrt(sandwich::vcovCL(c2_lm,
                            cluster = ~ state)["treatedTRUE", "treatedTRUE"]),
    "c1_analytic" = sqrt(
      varDD(ntx = nperstate_cohort1[1], 
            nctrl = nperstate_cohort1[-1],
            Tpre = Tpre, Tpost = Tpost,
            rho = with(state_varcors, rho[state != "tx02"]), 
            phi_t = with(state_varcors, phi[state != "tx02"]),
            phi_s = with(state_varcors, psi[state != "tx02"]),
            sigma_s = with(state_varcors, 
                           sqrt(outcome_var[state != "tx02"])))),
    "c2_analytic" = sqrt(
      varDD(ntx = nperstate_cohort2[1], 
            nctrl = nperstate_cohort2[-1],
            Tpre = Tpre, Tpost = Tpost,
            rho = with(state_varcors, rho[state != "tx01"]), 
            phi_t = with(state_varcors, phi[state != "tx01"]),
            phi_s = with(state_varcors, psi[state != "tx01"]),
            sigma_s = with(state_varcors, 
                           sqrt(outcome_var[state != "tx01"]))))
  )
  
  overall_estimates <- data.frame(
    "truth" = truth,
    "naive_est" = weighted.mean(
      x = with(did_estimates, c(c1_lm, c2_lm)),
      w = with(did_ses, c(1/c1_lm^2, 1/c2_lm^2))),
    "idadj_est" = weighted.mean(
      x = with(did_estimates, c(c1_lm, c2_lm)),
      w = with(did_ses, c(1/c1_lm_idadj^2, 1/c2_lm_idadj^2))),
    "stadj_est" = weighted.mean(
      x = with(did_estimates, c(c1_lm, c2_lm)),
      w = with(did_ses, c(1/c1_lm_stadj^2, 1/c2_lm_stadj^2))),
    "analytic_est" = weighted.mean(
      x = with(did_estimates, c(c1_lm, c2_lm)),
      w = with(did_ses, c(1/c1_analytic^2, 1/c2_analytic^2))),
    "naive_se" = 1 / sqrt(with(did_ses, c1_lm^(-2) + c2_lm^(-2))),
    "idadj_se" = 1 / sqrt(with(did_ses, c1_lm_idadj^(-2) + c2_lm_idadj^(-2))),
    "stadj_se" = 1 / sqrt(with(did_ses, c1_lm_stadj^(-2) + c2_lm_stadj^(-2))),
    "analytic_se" = 1 / sqrt(with(did_ses, 
                                  c1_analytic^(-2) + c2_analytic^(-2)))
    ) |>
    transform("naive_lb" = naive_est - qnorm(.975) * naive_se,
              "naive_ub" = naive_est + qnorm(.975) * naive_se,
              "idadj_lb" = idadj_est - qnorm(.975) * idadj_se,
              "idadj_ub" = idadj_est + qnorm(.975) * idadj_se,
              "stadj_lb" = stadj_est - qnorm(.975) * stadj_se,
              "stadj_ub" = stadj_est + qnorm(.975) * stadj_se,
              "analytic_lb" = analytic_est - qnorm(.975) * analytic_se,
              "analytic_ub" = analytic_est + qnorm(.975) * analytic_se) |> 
    transform("naive_coverage" = truth >= naive_lb & truth <= naive_ub,
              "idadj_coverage" = truth >= idadj_lb & truth <= idadj_ub,
              "stadj_coverage" = truth >= stadj_lb & truth <= stadj_ub,
              "analytic_coverage" = truth >= analytic_lb & truth <= analytic_ub)
  
  return(list("component_sums"  = component_sums,
              "component_sums_old" = component_sums_old,
              "component_means" = component_means,
              "did_estimates"   = did_estimates,
              "did_ses" = did_ses,
              "overall_estimates" = overall_estimates))
}
