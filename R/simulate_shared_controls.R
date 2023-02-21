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
#'   nctrlstates.
#' @param phi_t Within-period correlation, i.e., correlation between two
#'   observations from different people in the same state at the same time. If a
#'   single number, assumed constant across all states; else, a vector of length
#'   2 + nctrlstates.
#' @param phi_s Within-state correlation, i.e., correlation between two
#'   observations from different people in the same state at different times. If
#'   a single number, assumed constant across all states; else, a vector of
#'   length 2 + nctrlstates.
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
                                     phi_t = rep(0, 2 + nctrlstates), ## CHECK LENGTH
                                     phi_s = rep(0, 2 + nctrlstates), ## CHECK LENGTH
                                     tx_intshift_cohort1,
                                     tx_intshift_cohort2 = tx_intshift_cohort1,
                                     secular_cohort1 = function(x) 0,
                                     secular_cohort2 = secular_cohort1,
                                     state_fe_mean = 0,
                                     error_sd = rep(1, 2 + nctrlstates),
                                     verbose = FALSE) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(rho, add = checks)
  checkmate::assert_numeric(phi_t, add = checks)
  checkmate::assert_numeric(phi_s, add = checks)
  checkmate::reportAssertions(checks)
  
  if (length(nperstate_cohort1) == 1) {
    nperstate_cohort1 <- rep_len(nperstate_cohort1, nctrlstates + 1)
  } else if (length(nperstate_cohort1) != nctrlstates + 1)
    stop("nperstate_cohort1 must have length 1 or (nctrlstates + 1)")
  
  if (length(nperstate_cohort2) == 1) {
    nperstate_cohort2 <- rep_len(nperstate_cohort2, nctrlstates + 1)
  } else if (length(nperstate_cohort2) != nctrlstates + 1)
    stop("nperstate_cohort2 must have length 1 or (nctrlstates + 1)")
  
  if (length(nshared) == 1)
    nshared <- rep(nshared, nctrlstates)
  else if (length(nshared) != nctrlstates) {
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
                        txstatename = "tx1",
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
                        txstatename = "tx2",
                        cohort_name = "cohort2",
                        start_time = Delta)
  
  # Build matrix of variances for random effects
  state_varcors <- build_varcor_matrix(c1, c2, rho, phi_t, phi_s, error_sd)
  
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
  c2$cohort$id[c2$cohort$id %in% killed_ids] <- rep(shared_ids, each = 2)
  
  # Generate ID-specific random effects
  random_id_effect <- subset(rbind(c1$cohort, c2$cohort), 
                             select = c("id", "state"))
  random_id_effect <- random_id_effect[!duplicated(random_id_effect), ]
  random_id_effect <- 
    do.call(rbind, 
            lapply(split.data.frame(random_id_effect, random_id_effect$state), 
                   \(x) {
                     st <- as.character(unique(x$state))
                     x$rdmint_id <- rnorm(nrow(x), mean = 0,
                                          sd = with(state_varcors, 
                                                    sqrt(s2_id[state == st])))
                     x
                   }))
  
  random_time_effect <- rbind(expand.grid("state" = unique(c1$cohort$state),
                                          "time" = unique(c1$cohort$time)),
                              expand.grid("state" = unique(c2$cohort$state), 
                                          "time" = unique(c2$cohort$time)))
  random_time_effect <- random_time_effect[!duplicated(random_time_effect), ]
  random_time_effect$rdmint_time <- sapply(1:nrow(random_time_effect), \(i) {
    rnorm(n = 1, mean = 0, 
          sd = sqrt(state_varcors$s2_time[state_varcors$state == 
                                            random_time_effect$state[i]]))
  })
  
  data.frame("time" = sort(union(c1$cohort$time, c2$cohort$time)),
             "state" = union(c1$cohort$state, c2$cohort$state),
             "rdmint_time" = rnorm(Tpre + Tpost + Delta,
                                   mean = 0,
                                   sd = sqrt(state_varcors$s2_time[1])))
  
  state_fes <- expand.grid("state" = union(c1$cohort$state, c2$cohort$state),
                           # "cohort" = c(1, 2), 
                           stringsAsFactors = F)
  state_fes$state_eff <- MASS::mvrnorm(n = nrow(state_fes), mu = state_fe_mean, 
                                       Sigma = unique(state_varcors$s2_state),
                                       empirical = T)
  
  outcome_error <- 
    expand.grid(id = union(c1$cohort$id, c2$cohort$id),
                time = union(c1$cohort$time, c2$cohort$time))
  
  outcome_error$error <- rnorm(n = nrow(outcome_error),
                               mean = 0, sd = error_sd)
  
  # Merge ID-specific random effects and outcome error into cohorts
  c1$cohort <- merge(c1$cohort, random_id_effect, all.x = T)
  c2$cohort <- merge(c2$cohort, random_id_effect, all.x = T)
  c1$cohort <- merge(c1$cohort, random_time_effect, all.x = T)
  c2$cohort <- merge(c2$cohort, random_time_effect, all.x = T)
  c1$cohort <- merge(c1$cohort, state_fes, all.x = T)
  c2$cohort <- merge(c2$cohort, state_fes, all.x = T)
  
  c1$cohort <- merge(c1$cohort, outcome_error, all.x = T, by = c("id", "time"))
  c2$cohort <- merge(c2$cohort, outcome_error, all.x = T, by = c("id", "time"))
  
  # Merge in state fixed effects for early cohort (late cohort FEs depend on
  # early cohort!)
  c1$cohort$state_eff <- state_fes$state_eff[
    match(c1$cohort$state, state_fes$state)]
  
  # Merge in state fixed effects for late cohort (These values will change
  # depending on number of shared control individuals: shared folks will have a
  # different intercept than the state FE so their trajectories will stay the
  # same and not be cohort-dependent.)
  c2$cohort$state_eff <- state_fes$state_eff[
    match(c2$cohort$state, state_fes$state)]
  
  # Construct outcome variables in early cohort
  c1$cohort$Ymean <- with(c1$cohort,  
                          state_eff + 
                            secular_cohort1(time) +
                            tx_intshift_cohort1 * treated)
  
  c1$cohort$Y <- with(c1$cohort, Ymean + rdmint_id + rdmint_time + error)
  
  # # cohort 2 FEs if there is time overlap
  # if (Delta <= Tpre + Tpost) {
  #   adj_intercepts  <-
  #     aggregate(Ymean ~ state,  subset(c1$cohort, grepl("ctrl", state) &
  #                                        id %in% shared_ids & time == Delta),
  #               unique)
  #   adj_intercepts$disjoint <- 
  #     (nperstate_cohort2[-1] * 
  #        state_fes$state_eff[grepl("ctrl", state_fes$state) &
  #                              state_fes$cohort == 2] - 
  #        nshared * adj_intercepts$Ymean) / 
  #     (nperstate_cohort2[-1] - nshared)
  #   
  #   # Adjust late cohort's "fixed effects" so that we achieve the desired state
  #   # FE (the adjusted values will ensure the average is correct)
  #   c2$cohort <- 
  #     do.call(rbind, 
  #             c(lapply(split.data.frame(c2$cohort, ~ state), \(x) {
  #               st <- as.character(unique(x$state))
  #               if (grepl("tx", st)) {
  #                 x
  #               } else {
  #                 x$state_eff[x$id %in% shared_ids] <-
  #                   adj_intercepts$Ymean[adj_intercepts$state == st]
  #                 x$state_eff[!(x$id %in% shared_ids)] <-
  #                   adj_intercepts$disjoint[adj_intercepts$state == st]
  #                 x
  #               }
  #             }), make.row.names = F)
  #     )
  # }
  
  #CHECK CORRELATION
  # resids1 <- reshape(subset(c1$cohort, select = c("id", "time", "state", "Y", "Ymean")),
  #                    direction = "wide", idvar = c("id", "state"), timevar = "time")
  # resids1 <- transform(resids1, "residprod" = (Y.0 - Ymean.0) * (Y.1 - Ymean.1))
  # if (verbose)
  #   cat(paste("Within-person estimate: ", round(mean(resids1$residprod), 6), "\n"))
  # resids2 <- c1$cohort |> transform("resid" = Y - Ymean) |>
  #   subset(select = c("id", "time", "resid")) |>
  #    reshape(direction = "wide", idvar = "time", timevar = "id",
  #            v.names = "resid")
  #  within_period <- sum(apply(resids2[, sample(2:(ncol(resids2) - 1), 20000, replace = F)], 1, \(x) sum(combn(x[-1], 2, prod)))) /
  #   (nrow(resids2) * choose(ncol(resids2) - 1, 2))
  # if (verbose)
  #   cat(paste("Within-period estimate:", round(within_period, 6), "\n"))
  # resids3 <- t(resids2)
  # between_period <- sum(sapply(1:(ncol(resids3) - 1), \(j) {
  #   sum(sapply(1:nrow(resids3), \(i) {
  #     sum(resids3[i, j] * resids3[-i, (j + 1):ncol(resids3)])
  #   }))
  # })) / (nrow(resids3)
  
  c2$cohort$Ymean <- with(c2$cohort,
                          state_eff +
                            secular_cohort2(time) +
                            tx_intshift_cohort2 * treated)
  c2$cohort$Y <- with(c2$cohort, Ymean + rdmint_id + rdmint_time + error)
  
  # state_time_means <- structure(
  #   aggregate(Y ~ state + time, data = c1$cohort, mean),
  #   "names" = c("state", "time", "empMeanY")
  # )
  # c1$cohort <- merge(c1$cohort, state_time_means, by = c("state", "time"), 
  #                    all.y = T)
  # idcombn.cohort1.ctrl1 <- with(c1$cohort, t(combn(unique(id[state == "ctrl1"]), 2)))
  # test2 <- Reduce("+", lapply(idcombn.cohort1.ctrl1, \(x) {
  #   with(c1$cohort, crossprod(Y[id == x[1]] - empMeanY[id == x[1]],
  #                              Y[id == x[2]] - empMeanY[id == x[2]]))
  # }))
  
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
  
  c1disj   <- subset(c1$cohort, !(id %in% shared_ids) & state != "tx1")
  c1shared <- subset(c1$cohort, id %in% shared_ids & state != "tx1")
  c2disj   <- subset(c2$cohort, !(id %in% shared_ids) & state != "tx2")
  c2shared <- subset(c2$cohort, id %in% shared_ids & state != "tx2")
  
  # Linear Models
  c1_lm <- lm(Y ~ -1 + state + as.factor(time) + treated, data = c1$cohort)
  c2_lm <- lm(Y ~ -1 + state + as.factor(time) + treated, data = c2$cohort)
  
  # True overall treatment effect
  truth <- mean(c(tx_intshift_cohort1, tx_intshift_cohort2))
  
  component_sums <- sapply(c(c1_expand, c2_expand), \(x) sum(x$Y, na.rm = T))
  
  component_sums_old <- data.frame(
    "c1_tx_pre"    = with(c1$cohort, sum(Y[state == "tx1" & !post])),
    "c1_tx_post"   = with(c1$cohort, sum(Y[state == "tx1" & post])),
    "c1_disj_pre"  = with(c1disj,    sum(Y[!post])),
    "c1_disj_post" = with(c1disj,    sum(Y[post])),
    "c1_pre_dot"   = with(c1shared,  sum(Y[time %in% c1c2t$pre_dot])),
    "c1_pre_pre"   = with(c1shared,  sum(Y[time %in% c1c2t$pre_pre])),
    "c1_post_pre"  = with(c1shared,  sum(Y[time %in% c1c2t$post_pre])),
    "c1_post_post" = with(c1shared,  sum(Y[time %in% c1c2t$post_post])),
    "c1_post_dot"  = with(c1shared,  sum(Y[time %in% c1c2t$post_dot])),
    "c2_tx_pre"    = with(c2$cohort, sum(Y[state == "tx2" & !post])),
    "c2_tx_post"   = with(c2$cohort, sum(Y[state == "tx2" & post])),
    "c2_disj_pre"  = with(c2disj,    sum(Y[!post])),
    "c2_disj_post" = with(c2disj,    sum(Y[post])),
    "c2_pre_pre"   = with(c2shared,  sum(Y[time %in% c1c2t$pre_pre])),
    "c2_post_pre"  = with(c2shared,  sum(Y[time %in% c1c2t$post_pre])),
    "c2_post_post" = with(c2shared,  sum(Y[time %in% c1c2t$post_post])),
    "c2_dot_pre"   = with(c2shared,  sum(Y[time %in% c1c2t$dot_pre])),
    "c2_dot_post"  = with(c2shared,  sum(Y[time %in% c1c2t$dot_post]))
  )
  
  component_means <- data.frame(
    "c1_tx_pre"    = with(c1$cohort, mean(Y[state == "tx1" & !post])),
    "c1_tx_post"   = with(c1$cohort, mean(Y[state == "tx1" & post])),
    "c1_ctrl_pre"  = with(c1$cohort, mean(Y[state != "tx1" & !post])),
    "c1_ctrl_post" = with(c1$cohort, mean(Y[state != "tx1" & post])),
    "c2_tx_pre"    = with(c2$cohort, mean(Y[state == "tx2" & !post])),
    "c2_tx_post"   = with(c2$cohort, mean(Y[state == "tx2" & post])),
    "c2_ctrl_pre"  = with(c2$cohort, mean(Y[state != "tx2" & !post])),
    "c2_ctrl_post" = with(c2$cohort, mean(Y[state != "tx2" & post]))
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
    "c1_lm_adj" = 
      sqrt(sandwich::vcovCL(c1_lm, cluster = ~ id)["treatedTRUE", "treatedTRUE"]),
    "c2_lm_adj" = 
      sqrt(sandwich::vcovCL(c2_lm, cluster = ~ id)["treatedTRUE", "treatedTRUE"]),
    "c1_analytic" = sqrt(
      var_DD(Ntx = nperstate_cohort1[1], 
             nctrl = nperstate_cohort1[-1],
             Tpre = Tpre, Tpost = Tpost,
             rho = with(state_varcors, rho[state != "tx2"]), 
             phi_t = with(state_varcors, phi_t[state != "tx2"]),
             phi_s = with(state_varcors, phi_s[state != "tx2"]),
             sigma_s = with(state_varcors, 
                            sqrt(outcome_var[state != "tx2"])))),
    "c2_analytic" = sqrt(
      var_DD(Ntx = nperstate_cohort2[1], 
             nctrl = nperstate_cohort2[-1],
             Tpre = Tpre, Tpost = Tpost,
             rho = with(state_varcors, rho[state != "tx1"]), 
             phi_t = with(state_varcors, phi_t[state != "tx1"]),
             phi_s = with(state_varcors, phi_s[state != "tx1"]),
             sigma_s = with(state_varcors, 
                            sqrt(outcome_var[state != "tx1"]))))
  )
  
  overall_estimates <- data.frame(
    "truth" = truth,
    "est" = weighted.mean(x = with(did_estimates, c(c1_lm, c2_lm)),
                          w = with(did_ses, c(1/c1_lm^2, 1/c2_lm^2))),
    "naive_se" = 1 / sqrt(with(did_ses, c1_lm^(-2) + c2_lm^(-2))),
    "adj_se" = 1 / sqrt(with(did_ses, c1_lm_adj^(-2) + c2_lm_adj^(-2)))) |>
    transform("naive_lb" = est - 1.96 * naive_se,
              "naive_ub" = est + 1.96 * naive_se,
              "adj_lb" = est - 1.96 * adj_se,
              "adj_ub" = est + 1.96 * adj_se) |> 
    transform("naive_coverage" = truth >= naive_lb & truth <= naive_ub,
              "adj_coverage" = truth >= adj_lb & truth <= adj_ub) 
  
  return(list("component_sums"  = component_sums,
              "component_sums_old" = component_sums_old,
              "component_means" = component_means,
              "did_estimates"   = did_estimates,
              "did_ses" = did_ses,
              "overall_estimates" = overall_estimates))
}
