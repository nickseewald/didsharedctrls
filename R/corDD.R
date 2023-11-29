# corDD.R 
# Copyright 2023 Nicholas J. Seewald

# Contains user-accessible functions to compute correlation between
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


#' Estimate Correlation Due to Shared Control Individuals between
#' Difference-in-Differences Estimates
#'
#' @description In stacked difference-in-differences (DiD) analyses using
#'   individual-level data where non-policy-implementing units are used as
#'   comparators for multiple policy-implementing units, data from untreated
#'   individuals may be used across multiple analyses, thereby inducing
#'   correlation between effect estimates. This function, given information
#'   about the stacked DiD design, estimates this correlation between estimates.
#'
#' @param Ntx_cohort1 Number of individuals in the treated state in cohort 1
#' @param Ntx_cohort2 Number of individuals in the treated state in cohort 2
#' @param Ndisj_cohort1 Vector of numbers of unshared individuals in each of the
#'   control states who contribute to cohort 1. Length of this vector is taken
#'   to be the number of control states.
#' @param Ndisj_cohort2 Vector of numbers of unshared individuals in each of the
#'   control states who contribute to cohort 2. Must have the same length as
#'   `Ndisj_cohort1`.
#' @param Nshared Vector of numbers of shared control individuals in each of the
#'   control states who contribute to both cohorts 1 and 2. Must have the same
#'   length as `Ndisj_cohort1`.
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed the same across cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed the same across cohorts)
#' @param Delta Number of measurement occasions between cohort treatment times
#' @param rho Vector of within-person correlations for each state in the
#'   analysis. First element is assumed to be for the treated state. Recycled if
#'   length-1; otherwise must be length of `Ndisj_cohort1` + 1.
#' @param phi Vector of within-period correlations for each state in the
#'   analysis. First element is assumed to be for the treated state. Recycled if
#'   length-1; otherwise must be length of `Ndisj_cohort1` + 1. Must be
#'   elementwise less than rho.
#' @param psi Vector of within-state correlations for each state in the
#'   analysis. First element is assumed to be for the treated state. Recycled if
#'   length-1; otherwise must be length of `Ndisj_cohort1` + 1. Must be
#'   elementwise less than phi.
#' @param outcomeSD Vector of standard deviations of the outcome for each state
#'   in the analysis. First element is assumed to be for the treated state in
#'   cohort 1; second element the treated state in cohort 2.
#'
#' @return A scalar estimate of the correlation between the two cohorts' effect
#'   estimates.
#' @export
#'
#' @examples
#' corDD(Ntx_cohort1 = 500, Ntx_cohort2 = 450,
#'       Ndisj_cohort1 = c(250, 375, 400),
#'       Ndisj_cohort2 = c(175, 400, 350),
#'       Nshared = c(100, 200, 300),
#'       Tpre = 5, Tpost = 5, Delta = 3,
#'       rho = c(0.3, 0.29, 0.34, 0.22),
#'       phi = c(0.1, 0.07, 0.15, 0.11),
#'       psi = c(0.05, 0.02, 0.09, 0.05),
#'       outcomeSD = 2.5)
corDD <- function(Ntx_cohort1, Ntx_cohort2, Ndisj_cohort1, Ndisj_cohort2,
                   Nshared, Tpre, Tpost, Delta, rho, phi, psi, outcomeSD) {
  #### Check inputs
  checks <- checkmate::makeAssertCollection()
  # Only allow one treated state per cohort
  checkmate::assert_number(Ntx_cohort1, lower = 1, finite = T, add = checks)
  checkmate::assert_number(Ntx_cohort2, lower = 1, finite = T, add = checks)
  # Make sure Ndisj_cohort1 is a vector, the length of which determines
  # nctrlstates
  checkmate::assert_vector(Ndisj_cohort1, add = checks)
  nctrlstates <- length(Ndisj_cohort1)
  checkmate::assert_vector(Ndisj_cohort2, len = nctrlstates, add = checks)
  checkmate::assert_vector(Nshared,       len = nctrlstates, add = checks)
  
  # Recycle rho, phi, and psi
  # Order taken to be (tx state 1, tx state 2, ctrl state 1, ctrl state 2, ...)
  rho       <- expand_vec(rho, nctrlstates + 2)
  phi       <- expand_vec(phi, nctrlstates + 2)
  psi       <- expand_vec(psi, nctrlstates + 2)
  outcomeSD <- expand_vec(outcomeSD, nctrlstates + 2)
  
  checkmate::reportAssertions(checks)
  
  covDD(Ndisj_cohort1, Ndisj_cohort2, Nshared, Tpre, Tpost, Delta,
         rho[-(1:2)], phi[-(1:2)], psi[-(1:2)], outcomeSD[-(1:2)]) / 
    sqrt(varDD(Ntx_cohort1, nctrl = Ndisj_cohort1 + Nshared, Tpre, Tpost, 
                rho[-2], phi[-2], psi[-2], outcomeSD[-2]) *
           varDD(Ntx_cohort2, nctrl = Ndisj_cohort2 + Nshared, Tpre, Tpost,
                  rho[-1], phi[-1], psi[-1], outcomeSD[-1]))
}