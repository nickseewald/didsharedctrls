# auxiliary.R
# Copyright 2023 Nicholas J. Seewald

# Contains helper functions didsharedctrls

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


#' Compute Variances of Random Effects for Simulating Shared Controls Data
#'
#'
#' @param rho Exchangeable within-person correlation. If a single number,
#'   assumed constant across all states; else, a vector of length 2 + the number
#'   of control states, given in the order (cohort 1 treated, cohort 2 treated,
#'   control states).
#' @param phi Within-period correlation, i.e., correlation between two
#'   observations from different people in the same state at the same time. If a
#'   single number, assumed constant across all states; else, a vector of length
#'   2 + the number of control states, given in the order (cohort 1 treated,
#'   cohort 2 treated, control states).
#' @param psi Within-state correlation, i.e., correlation between two
#'   observations from different people in the same state at different times. If
#'   a single number, assumed constant across all states; else, a vector of
#'   length 2 + the number of control states, given in the order (cohort 1
#'   treated, cohort 2 treated, control states).
#' @param error_sd Standard deviation of additive error in data generative
#'   model. See Details.
#' @param c1 Optional cohort dataset produced by [generate_cohort()] used to
#'   extract state names in cohort 1 for export.
#' @param c2 Optional cohort dataset produced by [generate_cohort()] used to
#'   extract state names in cohort 2 for export.
#'
#' @return A `data.frame`. One row per state in the analysis, returning input
#'   correlation parameters as well as variances of random effects needed to
#'   produce the correct correlation structure for the between-estimate
#'   correlation formula in Seewald et al.
#'
#' @examples compute_raneff_vars(rho = c(.5, .2, .1))
compute_raneff_vars <- function(rho, phi, psi, error_sd, c1 = NULL, c2 = NULL) {
  if (is.null(c1) & is.null(c2)) {
    d <- data.frame("state" = NA)
  } else {
    d <- data.frame("state" = union(c1$cohort$state, c2$cohort$state))
  }
  
  # Add rho, phi, and psi to d
  d <- manage_cor_input(d, rho, "rho") |> 
    manage_cor_input(phi, "phi") |> 
    manage_cor_input(psi, "psi") |> 
    manage_cor_input(error_sd, "error_sd") |> 
    # Compute r, the ratio of psi to phi
    transform("r" = ifelse(psi == phi, 1, psi / phi)) |>
    # Compute individual and state-time-level random effect variances
    transform(
      "s2_id" = error_sd^2 * (rho - psi) / ((1 - rho) - (phi - psi)),
      "s2_sttime" = error_sd^2 * phi / ((1 - rho) - (phi - psi))
    ) |> 
    # Compute outcome variance
    transform(
      "outcome_var" = s2_id + s2_sttime + error_sd^2
    )
  
  # Throw error if any random effect variances are invalid given choices of rho,
  # phi, psi, error_sd.
  if (any(d[, c("s2_id", "s2_sttime")] < 0)) {
    dvec <- as.vector(d)
    str1 <- ifelse(any(dvec$s2_id < 0), 
                   paste("s2_id:", 
                         paste(dvec$state[dvec$s2_id < 0], collapse = ", "),
                         "\n"),
                   "")
    str2 <- ifelse(any(dvec$s2_sttime < 0), 
                   paste("s2_sttime:", 
                         paste(dvec$state[dvec$s2_sttime < 0], collapse = ", "),
                         "\n"),
                   "")
    stop(paste("Choices of rho, phi, and psi yield negative random",
               "intercept variances: ", str1, str2))
  }
  
  d
}


#' Construct a Correlation Matrix
#'
#' Given a parameter `rho`, build a `corstr`-type correlation matrix of
#' dimension `p`.
#'
#' @param rho Correlation parameter
#' @param p Dimension of desired matrix (will be \eqn{p * p})
#' @param corstr Correlation structure; one of "independence", "exchangeable",
#'   or "ar1". Partial matching used.
#' @return A \eqn{p * p} matrix of correlations.
#' @importFrom checkmate assert_numeric
#' @export
#'
#' @examples
#' corstr(rho = 0.3, p = 3, corstr = "exchangeable")
cormat <-
  function(rho,
           p,
           corstr = c("independence", "exchangeable", "ar1")) { 
    
    rho <- as.numeric(rho)
    checkmate::assert_numeric(rho, lower = -1, upper = 1)
    
    corstr <- match.arg(corstr)
    
    if (corstr == "independence") {
      diag(rep(1, p))
    } else if (corstr == "exchangeable") {
      m <- matrix(rep(rho, p ^ 2), nrow = p)
      diag(m) <- rep(1, p)
      m
    } else if (corstr == "ar1") {
      m <- diag(p)
      rho ^ (abs(row(m) - col(m)))
    } else if (corstr == "unstructured") {
      if (length(rho) != p)
        stop("for unstructured corstr, must have rho of length p")
      m <- diag(3)
      tpairs <- combn(1:p, 2)
      for (j in 1:ncol(tpairs)) {
        m[tpairs[1, j], tpairs[2, j]] <-
          m[tpairs[2, j], tpairs[1, j]] <- rho[j]
      }
      m
    }
  }

#' Recycle Scalars to Vectors of a Given Length
#'
#' Checks if a vector has the desired length and throws an error if not, or if
#' given a scalar argument, recycles it to the desired length.
#'
#' @param vec A vector or scalar
#' @param len Desired length
#'
#' @return A vector of length `len`.
#'
#' @examples
#' expand_vec(1, 5)
#' expand_vec(c(1, 2), 2)
#' expand_vec(c(1, 2), 3)
expand_vec <- function(vec, len) {
  if (length(vec) == 1) 
    vec <- rep(vec, len)
  checkmate::assert_vector(vec, len = len)
  vec
}

coverage_check <- function(est, se, df, param_value, conf.level = .95) {
  # Recycle param_value if needed
  param_value <- rep_len(param_value, length.out = length(est))
  tfactor <- qt(c((1 - conf.level) / 2, (1 + conf.level) / 2), df = df)
  ci <- est + se %o% tfactor
  sapply(1:length(param_value), \(i) param_value[i] > min(ci[i, ]) & 
           param_value[i] < max(ci[i, ]))
}

#' Convert Vector to `data.frame` for Generating Random Effect Variances
#'
#' Given a vector `r` of correlations (or other data), add a variable called
#' `rname` to a `data.frame` `d`.
#'
#' @param d A `data.frame` containing at least one variable called `state`
#'   indicating different units in a stacked DiD analysis (treated and control).
#' @param r Numeric vector to add to `d`
#' @param rname Name of variable to add to `d`
#'
#' @return A `data.frame`
#'
#' @examples x <- manage_cor_input(d = data.frame("state" = c("MD", "MI", "PA")), r = c(0.2, 0.3, 0.4), rname = "rho")
#' x <- manage_cor_input(d = x, r = c(0.1, 0.2, 0.3), rname = "phi")
#' x
manage_cor_input <- function(d, r, rname) {
  if (length(r) == 1) {
    d[[rname]] <- r
  } else if (length(r) != 1 & length(r) == nrow(d)) {
    if (is.null(names(r))) {
      d[[rname]] <- r
    } else if (length(setdiff(names(r), d$state)) == 0) {
      d[[rname]] <- r[match(names(r), d$state)]
    } else {
      stop(paste0("Names of ', ", rname, "' do not match generated state names.",
                  "Names must be", paste(d$state, collapse = ", ")))
    }
  } else {
    stop(paste0("'", rname, "' must be either length 1 or length ", nrow(d),
                " (the number of unique states)."))
  }
  d
}

#' Compute Durations of Time Overlap Intervals
#'
#' @description When two treated states' study periods overlap in time, the
#'   total duration of the time window covered by the two study periods (i.e.,
#'   from the start of the earlier period to the end of the later period) can be
#'   decomposed into a number of "overlap intervals" based on whether a study
#'   period has started, ended, or is pre- or post-treatment. This function
#'   computes the number of measurement occasions that fall into each of those
#'   possible overlap intervals.
#'
#' @param Tpre Number of measurement occasions in the pre-treatment period
#'   (assumed common to both cohorts)
#' @param Tpost Number of measurement occasions in the post-treatment period
#'   (assumed common to both cohorts)
#' @param Delta Number of measurement occasions between cohort 1's treatment 
#' time and cohort 2's
#' 
#' @importFrom zoo as.yearmon
#' @seealso \link{time_windows}
#
#' @returns 
#' A named list of overlap interval durations. Names are:
#' 
#' * `pre_dot`: Overlap interval in which the early cohort is in its 
#' pre-treatment period and the later cohort's study period has not yet started.
#' * `pre_pre`: Both cohorts are in their pre-treatment periods.
#' * `post_dot`: Early cohort is in its post-treatment period, later cohort's 
#' study period has not yet started.
#' * `post_pre`: Early cohort is in its post-treatment period, later cohort is 
#' in its pre-treatment period.
#' * `post_post`: Both cohorts are in their post-treatment periods.
#' * `dot_pre`: The early cohort's study period has ended, later cohort is in
#' its pre-treatment period.
#' * `dot_post`: The early cohort's study period has ended, later cohort is in
#' its post-treatment period.
#' 
#' Note that it is not possible for all overlap interval durations to be
#' non-zero: at least one overlap interval will not exist for any combination
#' of positive `Tpre`, `Tpost`, and `Delta`.
#' 
#' @export
#'
#' @examples
#' overlap_durations(Tpre = 3, Tpost = 3, Delta = 2)
overlap_durations <- function(Tpre, Tpost, Delta) {
  init1 <- zoo::as.yearmon("2010-01")
  pre1  <- seq(from = init1, by = 1/12, length.out = Tpre)
  tx1   <- max(pre1) + 1/12
  post1 <- seq(from = tx1, by = 1/12, length.out = Tpost)
  
  tx2 <- tx1 + Delta/12
  pre2 <- seq(to = tx2 - 1/12, by = 1/12, length.out = Tpre)
  post2 <- seq(from = tx2, by = 1/12, length.out = Tpost)
  
  pre_dot   <- zoo::as.yearmon(setdiff(pre1, union(pre2, post2)))
  pre_pre   <- zoo::as.yearmon(intersect(pre1, pre2))
  post_dot  <- zoo::as.yearmon(setdiff(post1, union(pre2, post2)))
  post_pre  <- zoo::as.yearmon(intersect(post1, pre2))
  post_post <- zoo::as.yearmon(intersect(post1, post2))
  dot_pre   <- zoo::as.yearmon(setdiff(pre2, union(pre1, post1)))
  dot_post  <- zoo::as.yearmon(setdiff(post2, union(pre1, post1)))
  
  return(list(pre_dot   = length(pre_dot),
              pre_pre   = length(pre_pre),
              post_dot  = length(post_dot),
              post_pre  = length(post_pre),
              post_post = length(post_post),
              dot_pre   = length(dot_pre),
              dot_post  = length(dot_post)))
}

#' Identify Stages of Overlapping Study Periods in Difference-in-Differences
#'
#' @param c1Times Vector of times in cohort 1's study period
#' @param c2Times Vector of times in cohort 2's study period
#' @param c1tStar First treated time in cohort 1. Must be an element of
#'   `c1Times`.
#' @param c2tStar First treated time in cohort 2. Must be an element of
#'   `c2Times`.
#'
#' @return A named list identifying which of `c1Times` and `c2Times` fall into
#'   each possible pair of pre/post periods.
#'   
#' @seealso \link{overlap_durations}
#'
#' @examples
#' time_windows(c1Times = c(0, 1, 2), c2Times = c(1, 2, 3), c1tStar = 1, c2tStar = 2)
time_windows <- function(c1Times, c2Times, c1tStar, c2tStar) {
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_subset(c1tStar, c1Times, empty.ok = FALSE, add = checks)
  checkmate::assert_subset(c2tStar, c2Times, empty.ok = FALSE, add = checks)
  checkmate::reportAssertions(checks)
  list(
    pre_dot   = c1Times[c1Times < min(c2Times) & c1Times < c1tStar],
    pre_pre   = c1Times[c1Times >= min(c2Times) & c1Times < c1tStar],
    post_dot  = c1Times[c1Times >= c1tStar & c1Times < min(c2Times)],
    post_pre  = c1Times[c1Times >= c1tStar & c1Times >= min(c2Times) & 
                          c1Times < min(c2tStar) & c1Times <= max(c1Times)],
    post_post = c1Times[c1Times >= c1tStar & c1Times >= c2tStar],
    dot_pre   = c2Times[c2Times > max(c1Times) & c2Times < c2tStar],
    dot_post  = c2Times[c2Times > max(c1Times) & c2Times >= c2tStar]
  )
}