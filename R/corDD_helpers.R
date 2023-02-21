# corDD_helpers.R 
# Copyright 2023 Nicholas J. Seewald

# Contains helper functions to compute correlation between
# difference-in-differences estimators due to shared control individuals

# This file is part of didstackcor
#
# didstackcor is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TKTK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with didstackcor If not, see <https://www.gnu.org/licenses/>.

#' Variance of Single-Treated-State Difference-in-Differences Estimator for the
#' Average ATT across Pre- and Post-Treatment Measurement Occasions
#'
#' @param ntx Number of individuals in the single treated state
#' @param nctrl Vector of numbers of individuals for each control state
#' @param Tpre Number of measurement occasions in the pre-treatment period
#' @param Tpost Number of measurement occasions in the post-treatment period
#' @param rho Vector of within-person correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param phi_t Vector of within-period correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param phi_s Vector of within-state correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param sigma_s Vector of standard deviations of the outcome for each state in
#'   the analysis. First element is assumed to be for the treated state in
#'   cohort 1; second element the treated state in cohort 2.
#'
#' @return A scalar
var_DD <- function(ntx, nctrl, Tpre, Tpost, rho, phi_t, phi_s, sigma_s) {
  checks <- checkmate::makeAssertCollection()
  
  checkmate::assert_vector(nctrl, add = checks)
  nctrlstates <- length(nctrl)
  
  rho     <- expand_vec(rho, nctrlstates + 1)
  phi_t   <- expand_vec(phi_t, nctrlstates + 1)
  phi_s   <- expand_vec(phi_s, nctrlstates + 1)
  sigma_s <- expand_vec(sigma_s, nctrlstates + 1)
  
  checkmate::reportAssertions(checks)
  
  var_tx_period(ntx, Tpre, rho[1], phi_t[1], phi_s[1], sigma_s[1]) +
    var_tx_period(ntx, Tpost, rho[1], phi_t[1], phi_s[1], sigma_s[1]) +
    var_ctrl_period(nctrl, Tpre, rho[-1], phi_t[-1], phi_s[-1], sigma_s[-1]) +
    var_ctrl_period(nctrl, Tpost, rho[-1], phi_t[-1], phi_s[-1], sigma_s[-1]) -
    2 * var_cov_tx_pre_post(ntx, Tpre, Tpost, rho[1], phi_s[1], sigma_s[1]) - 
    2 * var_cov_ctrl_pre_post(nctrl, Tpre, Tpost, rho[-1], phi_s[-1], sigma_s[-1])
}


#' Estimate Covariance Due to Shared Control Individuals between
#' Difference-in-Differences Estimates
#'
#' @param Ndisj_cohort1 Vector of numbers of unshared individuals in the control
#'   states in cohort 1
#' @param Ndisj_cohort2 Vector of numbers of unshared individuals in the control
#'   states in cohort 2
#' @param Nshared Vector of numbers of shared individuals in the control states
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed the same across cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed the same across cohorts)
#' @param Delta Number of measurement occasions between cohort treatment times
#' @param rho Vector of within-person correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param phi_t Vector of within-period correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param phi_s Vector of within-state correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param sigma_s Vector of standard deviations of the outcome for each state in
#'   the analysis. First element is assumed to be for the treated state in
#'   cohort 1; second element the treated state in cohort 2.
#'
#' @return A scalar estimate of the covariance between the two cohorts' effect
#'   estimates.
#'
#' @examples
#' 
covDD <- function(Ndisj_cohort1, Ndisj_cohort2, Nshared, Tpre, Tpost, Delta,
                   rho, phi_t, phi_s, sigma_s) {
  
  # Check inputs
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_vector(Ndisj_cohort1, add = checks)
  nctrlstates <- length(Ndisj_cohort1)
  checkmate::assert_vector(Ndisj_cohort2, len = nctrlstates, add = checks)
  checkmate::assert_vector(Nshared,       len = nctrlstates, add = checks)
  
  rho     <- expand_vec(rho,     nctrlstates)
  phi_t   <- expand_vec(phi_t,   nctrlstates)
  phi_s   <- expand_vec(phi_s,   nctrlstates)
  sigma_s <- expand_vec(sigma_s, nctrlstates)
  
  checkmate::reportAssertions(checks)
  
  Nctrl1 <- sum(Ndisj_cohort1 + Nshared)
  Nctrl2 <- sum(Ndisj_cohort2 + Nshared)
  
  # Compute covariance between post-period means
  post_post <- sum(sapply(1:nctrlstates, \(i) {
    post_post_cov(ndisj_cohort1 = Ndisj_cohort1[i], 
                  ndisj_cohort2 = Ndisj_cohort2[i],
                  nshared = Nshared[i], Tpre, Tpost, Delta,
                  rho = rho[i], phi_t = phi_t[i], phi_s = phi_s[i],
                  sigma_s = sigma_s[i])
  })) / (Nctrl1 * Nctrl2 * Tpost^2)
  
  # Compute covariance between pre-period means
  pre_pre <- sum(sapply(1:nctrlstates, \(i) {
    pre_pre_cov(ndisj_cohort1 = Ndisj_cohort1[i], 
                ndisj_cohort2 = Ndisj_cohort2[i],
                nshared = Nshared[i],  Tpre, Tpost, Delta,
                rho = rho[i], phi_t = phi_t[i], phi_s = phi_s[i],
                sigma_s = sigma_s[i]) 
  })) / (Nctrl1 * Nctrl2 * Tpre^2)
  
  # Compute covariance between cohort 1's post-period mean and cohort2's
  # pre-period mean
  post_pre <- sum(sapply(1:nctrlstates, \(i) {
    post_pre_cov(ndisj_cohort1 = Ndisj_cohort1[i],
                 ndisj_cohort2 = Ndisj_cohort2[i],
                 nshared = Nshared[i],  Tpre, Tpost, Delta,
                 rho = rho[i], phi_t = phi_t[i], phi_s = phi_s[i],
                 sigma_s = sigma_s[i]) 
  })) / (Nctrl1 * Nctrl2 * Tpre * Tpost)
  
  # Compute covariance between cohort 1's pre-period mean and cohort2's
  # post-period mean
  pre_post <- sum(sapply(1:nctrlstates, \(i) {
    pre_post_cov(ndisj_cohort1 = Ndisj_cohort1[i],
                 ndisj_cohort2 = Ndisj_cohort2[i],
                 nshared = Nshared[i],  Tpre, Tpost, Delta,
                 rho = rho[i], phi_t = phi_t[i], phi_s = phi_s[i],
                 sigma_s = sigma_s[i]) 
  })) / (Nctrl1 * Nctrl2 * Tpre * Tpost)
  
  # Combine above quantities to get total covariance
  post_post + pre_pre - post_pre - pre_post
}

#' Estimate the Correlation Due to Shared Control Individuals between
#' Pre-Treatment Means in a Single Control Unit in Difference-and-Differences
#'
#' @description In stacked difference-in-differences, treatment effect estimates
#'   for different 
#'
#' @param ndisj1 Number of unshared control individuals in the unit who
#'   contribute to cohort 1
#' @param ndisj2 Number of unshared control individuals in the unit who
#'   contribute to cohort 2
#' @param nshared Number of shared control individuals in the unit who
#'   contribute to both cohorts
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed the same across cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed the same across cohorts)
#' @param Delta Number of measurement occasions between cohort treatment times
#' @param rho Within-person correlation for individuals in the unit.
#' @param phi_t Within-period correlation for the unit, i.e., correlation
#'   between outcomes in different individuals at the same time.
#' @param phi_s Within-state correlation for the unit, i.e., correlation between
#'   outcomes in different individuals at different times.
#' @param sigma_s Standard deviation of the outcome in the unit.
#'
#' @return A scalar estimate of the covariance between two pre-treatment means
#'   in a difference-in-differences context.
#' 
pre_pre_cov <- function(ndisj_cohort1, ndisj_cohort2, nshared, Tpre, Tpost,
                        Delta, rho, phi_t, phi_s, sigma_s) {
  duration <- overlap_durations(Tpre, Tpost, Delta)
  durationC1 <- c(duration[grep("pre_", names(duration))],
                  "pre_disj" = Tpre)
  durationC2 <- c(duration[grep("_pre", names(duration))],
                  "disj_pre" = Tpre)
  
  combos <- expand.grid("C1" = names(durationC1), 
                        "C2" = names(durationC2), "cov" = NA)
  
  # Compute covariances with disjoint sums
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "disj_pre"] <-
    ndisj_cohort1 * ndisj_cohort2 * (duration$pre_pre * phi_t + 
                                       (Tpre^2 - duration$pre_pre) * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "dot_pre"] <-
    ndisj_cohort1 * nshared * (Tpre * duration$dot_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "pre_pre"] <-
    ndisj_cohort1 * nshared * (duration$pre_pre * phi_t +
                                 (duration$pre_pre * (Tpre - 1) * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "post_pre"] <-
    ndisj_cohort1 * nshared * (duration$post_pre * Tpre * phi_s) *
    sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "disj_pre"] <-
    ndisj_cohort2 * nshared * (duration$pre_pre * phi_t +
                                 duration$pre_pre * (Tpre - 1) * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "disj_pre"] <- 
    ndisj_cohort2 * nshared * (duration$pre_dot * Tpre * phi_s) * sigma_s^2
  
  # Compute covariances between sums of shared individuals
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "pre_pre"] <-
    nshared * with(duration, (pre_dot * pre_pre * rho) +
                     ((nshared - 1) * pre_dot * pre_pre * phi_s)) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "post_pre"] <- 
    nshared * with(duration, pre_dot * post_pre * rho +
                     (nshared - 1) * pre_dot * post_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "pre_pre"] <- 
    nshared * (duration$pre_pre + 
                 duration$pre_pre * (duration$pre_pre - 1) * rho +
                 (nshared - 1) * (duration$pre_pre * phi_t +
                                    duration$pre_pre * (duration$pre_pre - 1) *
                                    phi_s)) * 
    sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "dot_pre"] <-
    with(duration, nshared * 
           (pre_pre * dot_pre * rho +
              (nshared - 1) * pre_pre * dot_pre * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "post_pre"] <-
    with(duration, nshared * 
           (pre_pre * post_pre * rho +
              (nshared - 1) * pre_pre * post_pre * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "dot_pre"] <- 
    with(duration, nshared *
           (pre_dot * dot_pre * rho +
              (nshared - 1) * pre_dot * dot_pre * phi_s)) *
    sigma_s^2
  
  return(sum(combos$cov))
}

#' Estimate the Covariance  Due to Shared Control Individuals between
#' Post-Treatment Means in a Single Control Unit in Difference-and-Differences
#'
#' @description In stacked difference-in-differences in which control units are
#'   common across the stacked analyses, it is conceivable that some individuals
#'   in those units contribute to multiple analyses. This function estimates the
#'   correlation b
#'
#' @param ndisj1 Number of unshared control individuals in the unit who
#'   contribute to cohort 1
#' @param ndisj2 Number of unshared control individuals in the unit who
#'   contribute to cohort 2
#' @param nshared Number of shared control individuals in the unit who
#'   contribute to both cohorts
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed the same across cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed the same across cohorts)
#' @param Delta Number of measurement occasions between cohort treatment times
#' @param rho Within-person correlation for individuals in the unit.
#' @param phi_t Within-period correlation for the unit, i.e., correlation
#'   between outcomes in different individuals at the same time.
#' @param phi_s Within-state correlation for the unit, i.e., correlation between
#'   outcomes in different individuals at different times.
#' @param sigma_s Standard deviation of the outcome in the unit.
#'
#' @return A scalar estimate of the covariance between two pre-treatment means
#'   in a difference-in-differences context.
#'   
post_post_cov <- function(ndisj_cohort1, ndisj_cohort2, nshared, Tpre, Tpost, 
                          Delta, rho, phi_t, phi_s, sigma_s) {
  duration <- overlap_durations(Tpre, Tpost, Delta)
  durationC1 <- c(duration[grep("post_", names(duration))],
                  "post_disj" = Tpost)
  durationC2 <- c(duration[grep("_post", names(duration))],
                  "disj_post" = Tpost)
  
  combos <- expand.grid("C1" = names(durationC1), 
                        "C2" = names(durationC2), "cov" = NA)
  
  # Compute covariances with disjoint sums
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "disj_post"] <-
    ndisj_cohort1 * ndisj_cohort2 * (duration$post_post * phi_t + 
                                       (Tpost^2 - duration$post_post) * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "dot_post"] <-
    ndisj_cohort1 * nshared * (Tpost * duration$dot_post * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "post_post"] <-
    ndisj_cohort1 * nshared * (duration$post_post * phi_t +
                                 (duration$post_post * (Tpost - 1) * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "disj_post"] <-
    ndisj_cohort2 * nshared * (duration$post_dot * Tpost * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "disj_post"] <- 
    ndisj_cohort2 * nshared * (duration$post_pre * Tpost * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_post" & combos$C2 == "disj_post"] <-
    ndisj_cohort2 * nshared * (duration$post_post * phi_t +
                                 duration$post_post * (Tpost - 1) * phi_s) * sigma_s^2
  
  # Compute covariances between sums of shared individuals
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "post_post"] <-
    (nshared * with(duration, post_pre * post_post * rho) +
       (nshared * (nshared - 1) * with(duration, post_pre * post_post) *
          phi_s)) * sigma_s^2
  
  combos$cov[combos$C1 == "post_post" & combos$C2 == "post_post"] <- 
    nshared * (duration$post_post + 
                 duration$post_post * (duration$post_post - 1) * rho +
                 (nshared - 1) * (duration$post_post * phi_t +
                                    duration$post_post * (duration$post_post - 1) * phi_s)) * 
    sigma_s^2
  
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "post_post"] <- 
    nshared * with(duration, post_dot * post_post * rho +
                     (nshared - 1) * post_dot * post_post * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "dot_post"] <-
    with(duration, nshared * (
      post_pre * dot_post * rho +
        (nshared - 1) * post_pre * dot_post * phi_s
    )) * sigma_s^2
  
  combos$cov[combos$C1 == "post_post" & combos$C2 == "dot_post"] <-
    with(duration, nshared * 
           (post_post * dot_post * rho +
              (nshared - 1) * post_post * dot_post * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "dot_post"] <- 
    with(duration, nshared *
           (post_dot * dot_post * rho +
              (nshared - 1) * post_dot * dot_post * phi_s)) *
    sigma_s^2
  
  return(sum(combos$cov))
}

#' Estimate the Correlation Due to Shared Control Individuals between
#' Post- and Pre-Treatment Means in a Single Control Unit in 
#' Difference-and-Differences
#'
#' @description In stacked difference-in-differences with multip
#'
#' @param ndisj_cohort1 Number of unshared control individuals in the unit who
#'   contribute to cohort 1
#' @param ndisj_cohort2 Number of unshared control individuals in the unit who
#'   contribute to cohort 2
#' @param nshared Number of shared control individuals in the unit who
#'   contribute to both cohorts
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed the same across cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed the same across cohorts)
#' @param Delta Number of measurement occasions between cohort treatment times
#' @param rho Within-person correlation for individuals in the unit.
#' @param phi_t Within-period correlation for the unit, i.e., correlation
#'   between outcomes in different individuals at the same time.
#' @param phi_s Within-state correlation for the unit, i.e., correlation between
#'   outcomes in different individuals at different times.
#' @param sigma_s Standard deviation of the outcome in the unit.
#'
#' @return A scalar estimate of the covariance between two pre-treatment means
#'   in a difference-in-differences context.
#'
post_pre_cov <- function(ndisj_cohort1, ndisj_cohort2, nshared, Tpre, Tpost,
                         Delta, rho, phi_t, phi_s, sigma_s) {
  duration <- overlap_durations(Tpre, Tpost, Delta)
  durationC1 <- c(duration[grep("post_", names(duration))],
                  "post_disj" = Tpost)
  durationC2 <- c(duration[grep("_pre", names(duration))],
                  "disj_pre" = Tpre)
  
  combos <- expand.grid("C1" = names(durationC1), 
                        "C2" = names(durationC2), "cov" = NA)
  
  # Compute covariances with disjoint sums
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "disj_pre"] <-
    ndisj_cohort1 * ndisj_cohort2 * (duration$post_pre * phi_t + 
                                       (Tpost * Tpre - duration$post_pre) * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "dot_pre"] <-
    ndisj_cohort1 * nshared * (Tpost * duration$dot_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "post_pre"] <-
    ndisj_cohort1 * nshared * (duration$post_pre * phi_t +
                                 (duration$post_pre * (Tpost - 1) * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "pre_pre"] <-
    ndisj_cohort1 * nshared * (Tpost * duration$pre_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "disj_pre"] <-
    ndisj_cohort2 * nshared * (duration$post_dot * Tpre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_post" & combos$C2 == "disj_pre"] <- 
    ndisj_cohort2 * nshared * (duration$post_post * Tpre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "disj_pre"] <-
    ndisj_cohort2 * nshared * (duration$post_pre * phi_t +
                                 duration$post_pre * (Tpre - 1) * phi_s) * sigma_s^2
  
  # Compute covariances between sums of shared individuals
  
  ## This combo cannot exist: only one of these can be non-zero duration
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "pre_pre"] <- 0
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "pre_pre"] <- 
    nshared * with(duration, post_pre * pre_pre * rho +
                     (nshared - 1) * post_pre * pre_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_post" & combos$C2 == "pre_pre"] <- 
    nshared * with(duration, post_post * pre_pre * rho +
                     (nshared - 1) * post_post * pre_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "post_pre"] <-
    nshared * with(duration, post_dot * post_pre * rho +
                     (nshared - 1) * post_dot * post_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "post_pre"] <-
    nshared * with(duration, post_pre + post_pre * (post_pre - 1) * rho +
                     (nshared - 1) * (post_pre * phi_t +
                                        post_pre * (post_pre - 1) * phi_s)) *
    sigma_s^2
  
  combos$cov[combos$C1 == "post_post" & combos$C2 == "post_pre"] <- 
    nshared * with(duration, post_post * post_pre * rho +
                     (nshared - 1) * post_post * post_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_dot" & combos$C2 == "dot_pre"] <- 
    nshared * with(duration, post_dot * dot_pre * rho +
                     (nshared - 1) * post_dot * dot_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "post_pre" & combos$C2 == "dot_pre"] <- 
    nshared * with(duration, post_pre * dot_pre * rho +
                     (nshared - 1) * post_pre * dot_pre * phi_s) * sigma_s^2
  
  ## This combo cannot exist: only one of these can be non-zero duration
  combos$cov[combos$C1 == "post_post" & combos$C2 == "dot_pre"] <- 0
  
  return(sum(combos$cov))
  # combos
}

#' Estimate the Correlation Due to Shared Control Individuals between
#' Pre- and Post-Treatment Means in a Single Control Unit in 
#' Difference-and-Differences
#'
#' @description In stacked difference-in-differences with multip
#'
#' @param ndisj1 Number of unshared control individuals in the unit who
#'   contribute to cohort 1
#' @param ndisj2 Number of unshared control individuals in the unit who
#'   contribute to cohort 2
#' @param nshared Number of shared control individuals in the unit who
#'   contribute to both cohorts
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed the same across cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed the same across cohorts)
#' @param Delta Number of measurement occasions between cohort treatment times
#' @param rho Within-person correlation for individuals in the unit.
#' @param phi_t Within-period correlation for the unit, i.e., correlation
#'   between outcomes in different individuals at the same time.
#' @param phi_s Within-state correlation for the unit, i.e., correlation between
#'   outcomes in different individuals at different times.
#' @param sigma_s Standard deviation of the outcome in the unit.
#'
#' @return A scalar estimate of the covariance between two pre-treatment means
#'   in a difference-in-differences context.
#'
pre_post_cov <- function(ndisj_cohort1, ndisj_cohort2, nshared, Tpre, Tpost,
                         Delta, rho, phi_t, phi_s, sigma_s) {
  duration <- overlap_durations(Tpre, Tpost, Delta)
  durationC1 <- c(duration[grep("pre_", names(duration))],
                  "pre_disj" = Tpre)
  durationC2 <- c(duration[grep("_post", names(duration))],
                  "disj_post" = Tpost)
  
  combos <- expand.grid("C1" = names(durationC1),
                        "C2" = names(durationC2), "cov" = NA)
  
  # Compute covariances with disjoint sums
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "disj_post"] <-
    ndisj_cohort1 * ndisj_cohort2 * Tpre * Tpost * phi_s * sigma_s^2
  
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "dot_post"] <-
    ndisj_cohort1 * nshared * Tpre * duration$dot_post * phi_s * sigma_s^2
  
  combos$cov[combos$C1 == "pre_disj" & combos$C2 == "post_post"] <-
    ndisj_cohort1 * nshared * Tpre * duration$post_post * phi_s * sigma_s^2
  
  combos$cov[combos$C1 == "post_disj" & combos$C2 == "pre_pre"] <-
    ndisj_cohort1 * nshared * (Tpost * duration$pre_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "disj_post"] <-
    ndisj_cohort2 * nshared * duration$pre_pre * Tpost * phi_s * sigma_s^2
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "disj_post"] <- 
    ndisj_cohort2 * nshared * duration$pre_dot * Tpost * phi_s * sigma_s^2
  
  # Compute covariances between sums of shared individuals
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "post_post"] <- 
    nshared * with(duration, pre_dot * post_post * rho +
                     (nshared - 1) * pre_dot * post_post * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "post_post"] <- 
    nshared * with(duration, post_post * pre_pre * rho +
                     (nshared - 1) * post_post * pre_pre * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_dot" & combos$C2 == "dot_post"] <-
    nshared * with(duration, pre_dot * dot_post * rho +
                     (nshared - 1) * pre_dot * dot_post * phi_s) * sigma_s^2
  
  combos$cov[combos$C1 == "pre_pre" & combos$C2 == "dot_post"] <-
    nshared * with(duration, pre_pre * dot_post * rho +
                     (nshared - 1) * (pre_pre * dot_post * phi_s)) *
    sigma_s^2
  
  return(sum(combos$cov))
}

#' Compute Variance of Mean Outcome in Treated State in either Pre or
#' Post-Treatment Period
#'
#' @param ntx Number of individuals in the treated unit
#' @param nTimes Number of measurement occasions in the period (pre or post tx)
#' @param rho Within-person correlation
#' @param phi_t Within-period correlation
#' @param phi_s Between-period correlation
#' @param sigma_tx Variance of outcome in treated unit
#'
#' @return A scalar
var_tx_period <- function(ntx, nTimes, rho, phi_t, phi_s, sigma_tx) {
  (ntx * (nTimes + nTimes * (nTimes - 1) * rho) +
     ntx * (ntx - 1) * (nTimes * phi_t + nTimes * (nTimes - 1) * phi_s)) *
    sigma_tx^2 / (ntx^2 * nTimes^2)
}

#' Compute Variance of Mean Outcome Averaged over Control States in either Pre
#' or Post-Treatment Period
#'
#' @param nctrl Vector of numbers of individuals in each control state
#' @param nTimes Number of measurement occasions in the period (pre or post tx)
#' @param rho Vector of within-person correlations for each control state
#' @param phi_t Vector of within-period correlations for each control state
#' @param phi_s Vector of between-period correlations for each control state
#' @param sigma_s Vector of outcome variances in each control state
#' 
#' @return
var_ctrl_period <- function(nctrl, nTimes, rho, phi_t, phi_s, sigma_s) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_vector(nctrl, add = checks)
  checkmate::assert_vector(rho,     len = length(nctrl), add = checks)
  checkmate::assert_vector(phi_t,   len = length(nctrl), add = checks)
  checkmate::assert_vector(phi_s,   len = length(nctrl), add = checks)
  checkmate::assert_vector(sigma_s, len = length(nctrl), add = checks)
  checkmate::reportAssertions(checks)
  
  sum(sapply(1:length(nctrl), \(i) {
    (nctrl[i] * (nTimes + nTimes * (nTimes - 1) * rho[i]) + 
       nctrl[i] * (nctrl[i] - 1) * 
       (nTimes * phi_t[i] + nTimes * (nTimes - 1) * phi_s[i])) * sigma_s[i]^2
  })) / (sum(nctrl)^2 * nTimes^2)
}

#' Compute Covariance between Pre- and Post-Treatment Mean Outcomes in a Treated State
#'
#' @param ntx Number of individuals in the treated unit
#' @param Tpre Number of measurement occasions in the pre-treatment period
#' @param Tpost Number of measurement occasions in the post-treatment period
#' @param rho Within-person correlation
#' @param phi_s Between-period correlation
#' @param sigma_tx Variance of outcome in treated unit
#'
#' @return A scalar
var_cov_tx_pre_post <- function(ntx, Tpre, Tpost, rho, phi_s, sigma_tx) {
  (ntx * Tpre * Tpost * rho + ntx * (ntx - 1) * Tpre * Tpost * phi_s) *
    sigma_tx^2 / (ntx^2 * Tpre * Tpost)
}

#' Compute Covariance between Pre- and Post-Treatment Mean Outcomes Averaged
#' over All Control States
#'
#' @param nctrl Vector of numbers of individuals in each control state
#' @param Tpre Number of measurement occasions in the pre-treatment period
#' @param Tpost Number of measurement occasions in the post-treatment period
#' @param rho Vector of within-person correlations for each control state
#' @param phi_s Vector of between-period correlations for each control state
#' @param sigma_s Vector of outcome variances in each control state
#'
#' @return
var_cov_ctrl_pre_post <- function(nctrl, Tpre, Tpost, rho, phi_s, sigma_s) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_vector(nctrl, add = checks)
  checkmate::assert_vector(rho,     len = length(nctrl), add = checks)
  checkmate::assert_vector(phi_s,   len = length(nctrl), add = checks)
  checkmate::assert_vector(sigma_s, len = length(nctrl), add = checks)
  checkmate::reportAssertions(checks)
  
  sum(sapply(1:length(nctrl), \(i) {
    (nctrl[i] * Tpre * Tpost * rho[i] + 
       nctrl[i] * (nctrl[i] - 1) * Tpre * Tpost * phi_s[i]) * sigma_s[i]^2
  })) / (sum(nctrl)^2 * Tpre * Tpost)
}