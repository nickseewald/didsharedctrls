# corDD.R 
# Copyright 2023 Nicholas J. Seewald

# Contains user-accessible functions to compute correlation between
# difference-in-differences estimators due to shared control individuals

# This file is part of didstackcor
#
# didstackcor is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# didstackcor is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with didstackcor If not, see <https://www.gnu.org/licenses/>.


#' Estimate Correlation Due to Shared Control Individuals between
#' Difference-in-Differences Estimates
#'
#' @description In stacked difference-in-differences in which control units are
#'   reused across
#'
#' @param Ntx_cohort1 Number of individuals in the treated state in cohort 1
#' @param Ntx_cohort2 Number of individuals in the treated state in cohort 2
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
#' @return A scalar estimate of the correlation between the two cohorts' effect
#'   estimates.
#' @export
#'
#' @examples
corDD <- function(Ntx_cohort1, Ntx_cohort2, Ndisj_cohort1, Ndisj_cohort2,
                   Nshared, Tpre, Tpost, Delta, rho, phi_t, phi_s, sigma_s) {
  # Check inputs
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_vector(Ndisj_cohort1, add = checks)
  nctrlstates <- length(Ndisj_cohort1)
  checkmate::assert_vector(Ndisj_cohort2, len = nctrlstates, add = checks)
  checkmate::assert_vector(Nshared,       len = nctrlstates, add = checks)
  
  rho     <- expand_vec(rho, nctrlstates + 1)
  phi_t   <- expand_vec(phi_t, nctrlstates + 1)
  phi_s   <- expand_vec(phi_s, nctrlstates + 1)
  sigma_s <- expand_vec(sigma_s, nctrlstates + 2)
  
  checkmate::reportAssertions(checks)
  
  covDD(Ndisj_cohort1, Ndisj_cohort2, Nshared, Tpre, Tpost, Delta,
         rho[-1], phi_t[-1], phi_s[-1], sigma_s[-(1:2)]) / 
    sqrt(var_DD(Ntx_cohort1, nctrl = Ndisj_cohort1 + Nshared, Tpre, Tpost, 
                rho, phi_t, phi_s, sigma_s[-2]) *
           var_DD(Ntx_cohort2, nctrl = Ndisj_cohort2 + Nshared, Tpre, Tpost,
                  rho, phi_t, phi_s, sigma_s[-1]))
}