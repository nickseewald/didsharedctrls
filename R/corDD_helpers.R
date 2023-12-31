# corDD_helpers.R 
# Copyright 2023 Nicholas J. Seewald

# Contains helper functions to compute correlation between
# difference-in-differences estimators due to shared control individuals

# This file is part of didsharedctrls
#
# didsharedctrls is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# didsharedctrls is distributed in the hope that it will be useful,
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
#' @param phi Vector of within-period correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param psi Vector of within-state correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param outcomeSD Vector of standard deviations of the outcome for each state in
#'   the analysis. First element is assumed to be for the treated state in
#'   cohort 1; second element the treated state in cohort 2.
#'
#' @return A scalar
varDD <- function(ntx, nctrl, Tpre, Tpost, rho, phi, psi, outcomeSD) {
  
  checkmate::assert_vector(nctrl)

  nctrlstates <- length(nctrl)
  
  rho     <- expand_vec(rho, nctrlstates + 1)
  phi     <- expand_vec(phi, nctrlstates + 1)
  psi     <- expand_vec(psi, nctrlstates + 1)
  outcomeSD <- expand_vec(outcomeSD, nctrlstates + 1)
  
  var_tx_period(ntx, Tpre, rho[1], phi[1], psi[1], outcomeSD[1]) +
    var_tx_period(ntx, Tpost, rho[1], phi[1], psi[1], outcomeSD[1]) +
    var_ctrl_period(nctrl, Tpre, rho[-1], phi[-1], psi[-1], outcomeSD[-1]) +
    var_ctrl_period(nctrl, Tpost, rho[-1], phi[-1], psi[-1], outcomeSD[-1]) -
    2 * var_cov_tx_pre_post(ntx, Tpre, Tpost, rho[1], psi[1], outcomeSD[1]) - 
    2 * var_cov_ctrl_pre_post(nctrl, Tpre, Tpost, rho[-1], psi[-1], outcomeSD[-1])
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
#' @param phi Vector of within-period correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param psi Vector of within-state correlations for each state in the
#'   analysis. First element is assumed to be for the treated state.
#' @param sigma_s Vector of standard deviations of the outcome for each state in
#'   the analysis. First element is assumed to be for the treated state in
#'   cohort 1; second element the treated state in cohort 2.
#'
#' @return A scalar estimate of the covariance between the two cohorts' effect
#'   estimates.
#' 
covDD <- function(Ndisj_cohort1, Ndisj_cohort2, Nshared, Tpre, Tpost, Delta,
                   rho, phi, psi, outcomeSD) {
  
  # Function is not exported, so will assume all input is well-formed
  
  # Get total number of control individuals per state per cohort
  Ncohort1 <- Ndisj_cohort1 + Nshared
  Ncohort2 <- Ndisj_cohort2 + Nshared
  
  timeFactor <- (Tpre * Tpost) ^ (-2) *
    (Tpre ^ 2 * max(Tpost - Delta, 0) +
       Tpost ^ 2 * max(Tpre - Delta, 0) -
       Tpre * Tpost * min(c(Tpre, Tpost, Delta, max(Tpre + Tpost - Delta, 0))))
  
  s <- ((Ncohort1 * (phi - psi)) * Ncohort2 + 
          Nshared * (1 - rho - (phi - psi))) * outcomeSD^2
  
  # s <- (Ncohort1 * Ncohort2 * (phi - psi) + 
  #         Nshared * (1 - rho - (phi - psi))) * sigma_s^2
  
  timeFactor * (sum(s) / sum(Ncohort1)) / sum(Ncohort2)
}

#' Compute Variance of Mean Outcome in Treated State in either Pre or
#' Post-Treatment Period
#'
#' @param ntx Number of individuals in the treated unit
#' @param nTimes Number of measurement occasions in the period (pre or post tx)
#' @param rho Within-person correlation
#' @param phi Within-period correlation
#' @param psi Between-period correlation
#' @param outcomeSD Variance of outcome in treated unit
#'
#' @return A scalar
var_tx_period <- function(ntx, nTimes, rho, phi, psi, outcomeSD) {
  (ntx * nTimes)^(-2) * (
    ntx * nTimes * (1 + (nTimes - 1) * rho) +
      ntx * (ntx - 1) * nTimes * (phi + (nTimes - 1) * psi)
  ) * outcomeSD^2
}

#' Compute Variance of Mean Outcome Averaged over Control States in either Pre-
#' or Post-Treatment Period
#'
#' @param nctrl Vector of numbers of individuals in each control state
#' @param nTimes Number of measurement occasions in the period (pre or post tx)
#' @param rho Vector of within-person correlations for each control state
#' @param phi Vector of within-period correlations for each control state
#' @param psi Vector of between-period correlations for each control state
#' @param outcomeSD Vector of outcome variances in each control state
#' 
#' @return A scalar, the variance of the mean outcome in control states over 
#' either the pre- or post-treatment period
var_ctrl_period <- function(nctrl, nTimes, rho, phi, psi, outcomeSD) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_vector(nctrl, add = checks)
  checkmate::assert_vector(rho,     len = length(nctrl), add = checks)
  checkmate::assert_vector(phi,     len = length(nctrl), add = checks)
  checkmate::assert_vector(psi,     len = length(nctrl), add = checks)
  checkmate::assert_vector(outcomeSD, len = length(nctrl), add = checks)
  checkmate::reportAssertions(checks)
  
  sum(sapply(1:length(nctrl), \(i) {
    (nctrl[i] * (nTimes + nTimes * (nTimes - 1) * rho[i]) + 
       nctrl[i] * (nctrl[i] - 1) * 
       (nTimes * phi[i] + nTimes * (nTimes - 1) * psi[i])) * outcomeSD[i]^2
  })) / (sum(nctrl)^2 * nTimes^2)
}

#' Compute Covariance between Pre- and Post-Treatment Mean Outcomes in a Treated State
#'
#' @param ntx Number of individuals in the treated unit
#' @param Tpre Number of measurement occasions in the pre-treatment period
#' @param Tpost Number of measurement occasions in the post-treatment period
#' @param rho Within-person correlation
#' @param psi Between-period correlation
#' @param outcomeSD Variance of outcome in treated unit
#'
#' @return A numeric scalar
var_cov_tx_pre_post <- function(ntx, Tpre, Tpost, rho, psi, outcomeSD) {
  (ntx * Tpre * Tpost * rho + ntx * (ntx - 1) * Tpre * Tpost * psi) *
    outcomeSD^2 / (ntx^2 * Tpre * Tpost)
}

#' Compute Covariance between Pre- and Post-Treatment Mean Outcomes Averaged
#' over All Control States
#'
#' @param nctrl Vector of numbers of individuals in each control state
#' @param Tpre Number of measurement occasions in the pre-treatment period
#' @param Tpost Number of measurement occasions in the post-treatment period
#' @param rho Vector of within-person correlations for each control state
#' @param psi Vector of between-period correlations for each control state
#' @param outcomeSD Vector of outcome variances in each control state
#'
#' @return A numeric scalar
var_cov_ctrl_pre_post <- function(nctrl, Tpre, Tpost, rho, psi, outcomeSD) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_vector(nctrl, add = checks)
  checkmate::assert_vector(rho,     len = length(nctrl), add = checks)
  checkmate::assert_vector(psi,     len = length(nctrl), add = checks)
  checkmate::assert_vector(outcomeSD, len = length(nctrl), add = checks)
  checkmate::reportAssertions(checks)
  
  sum(sapply(1:length(nctrl), \(i) {
    (nctrl[i] * Tpre * Tpost * rho[i] + 
       nctrl[i] * (nctrl[i] - 1) * Tpre * Tpost * psi[i]) * outcomeSD[i]^2
  })) / (sum(nctrl)^2 * Tpre * Tpost)
}