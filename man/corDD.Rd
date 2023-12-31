% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corDD.R
\name{corDD}
\alias{corDD}
\title{Estimate Correlation Due to Shared Control Individuals between
Difference-in-Differences Estimates}
\usage{
corDD(
  Ntx_cohort1,
  Ntx_cohort2,
  Ndisj_cohort1,
  Ndisj_cohort2,
  Nshared,
  Tpre,
  Tpost,
  Delta,
  rho,
  phi,
  psi,
  outcomeSD
)
}
\arguments{
\item{Ntx_cohort1}{Number of individuals in the treated state in cohort 1}

\item{Ntx_cohort2}{Number of individuals in the treated state in cohort 2}

\item{Ndisj_cohort1}{Vector of numbers of unshared individuals in each of the
control states who contribute to cohort 1. Length of this vector is taken
to be the number of control states.}

\item{Ndisj_cohort2}{Vector of numbers of unshared individuals in each of the
control states who contribute to cohort 2. Must have the same length as
\code{Ndisj_cohort1}.}

\item{Nshared}{Vector of numbers of shared control individuals in each of the
control states who contribute to both cohorts 1 and 2. Must have the same
length as \code{Ndisj_cohort1}.}

\item{Tpre}{Number of measurement occasions in the pre-treatment period
(assumed the same across cohorts)}

\item{Tpost}{Number of measurement occasions in the post-treatment period
(assumed the same across cohorts)}

\item{Delta}{Number of measurement occasions between cohort treatment times}

\item{rho}{Vector of within-person correlations for each state in the
analysis. First element is assumed to be for the treated state. Recycled if
length-1; otherwise must be length of \code{Ndisj_cohort1} + 1.}

\item{phi}{Vector of within-period correlations for each state in the
analysis. First element is assumed to be for the treated state. Recycled if
length-1; otherwise must be length of \code{Ndisj_cohort1} + 1. Must be
elementwise less than rho.}

\item{psi}{Vector of within-state correlations for each state in the
analysis. First element is assumed to be for the treated state. Recycled if
length-1; otherwise must be length of \code{Ndisj_cohort1} + 1. Must be
elementwise less than phi.}

\item{outcomeSD}{Vector of standard deviations of the outcome for each state
in the analysis. First element is assumed to be for the treated state in
cohort 1; second element the treated state in cohort 2.}
}
\value{
A scalar estimate of the correlation between the two cohorts' effect
estimates.
}
\description{
In stacked difference-in-differences (DiD) analyses using
individual-level data where non-policy-implementing units are used as
comparators for multiple policy-implementing units, data from untreated
individuals may be used across multiple analyses, thereby inducing
correlation between effect estimates. This function, given information
about the stacked DiD design, estimates this correlation between estimates.
}
\examples{
corDD(Ntx_cohort1 = 500, Ntx_cohort2 = 450,
      Ndisj_cohort1 = c(250, 375, 400),
      Ndisj_cohort2 = c(175, 400, 350),
      Nshared = c(100, 200, 300),
      Tpre = 5, Tpost = 5, Delta = 3,
      rho = c(0.3, 0.29, 0.34, 0.22),
      phi = c(0.1, 0.07, 0.15, 0.11),
      psi = c(0.05, 0.02, 0.09, 0.05),
      outcomeSD = 2.5)
}
