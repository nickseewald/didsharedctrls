# generate_cohort.R
# Copyright 2023 Nicholas J. Seewald, PhD

#' Generate a Data Skeleton for Staggered Adoption Difference in Differences
#' 
#' Create a "long" data frame to serve as a base for simulation studies for
#' staggered adoption difference in differences methods.
#'
#' @param nperstate number of individuals per state; can be a single integer or 
#' a vector of length `ntxstates` + `nctrlstates`. If a single integer, it is 
#' assumed that each state has `nperstate` individuals.
#' @param nctrlstates number of control states
#' @param Tobs number of observation periods, assumed to be in months
#' @param Tpost Number of observations in the study period for which the
#' treated state is treated. Must be less than `monthly_obs`.
#'
#' @return
#' @export
#'
#' @examples
generate_cohort <- function(nperstate, nctrlstates,
                            Tobs, Tpost,
                            txstatename = "tx1",
                            cohort_name = txstatename,
                            start_time = 0) {
  
  checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(nctrlstates, lower = 1, add = checks)
  checkmate::assert_integerish(nperstate, lower = 1, 
                               max.len = nctrlstates + 1, add = checks)
  checkmate::assert_integerish(Tobs, lower = 2, add = checks)
  checkmate::assert_integerish(Tpost, upper = Tobs - 1, add = checks)
  checkmate::assert_character(txstatename, add = checks)
  checkmate::reportAssertions(checks)
  
  cohortStates <- c(txstatename, 
                    paste0(c(rep("ctrl", nctrlstates)),
                           formatC(1:nctrlstates, width = nchar(nctrlstates),
                                   flag = "0")))
  
  if (length(nperstate) == 1) {
    nperstate <- rep(nperstate, 1 + nctrlstates)
  }
  
  names(nperstate) <- cohortStates
  
  
  crc <- digest::getVDigest(algo = "crc32")
  
  cohort <-
    do.call(rbind,
            lapply(
              cohortStates,
              \(st) {
                expand.grid(
                  "time" = start_time:(Tobs - 1 + start_time),
                  "personIndex" =
                    1:nperstate[which(cohortStates == st)],
                  "state" = st
                ) |> 
                  transform(
                    treated = as.logical(
                      grepl("tx", st) * (time >= Tobs + start_time - Tpost)),
                    txgroup = grepl("tx", st) * Tpost + 
                      (1 - grepl("tx", st)) * 999
                  )
              })) |>
    transform(cohort = cohort_name,
              month = time %% 12,
              year = time %/% 12,
              txgroup = ifelse(grepl("tx", state), 
                               Tpost[match(state, cohortStates)], 999),
              txtime = pmax(time - txgroup + 1, 0),
              id = crc(paste0(cohort_name, state, personIndex)),
              post = time >= start_time + (Tobs - Tpost)
    ) |> 
    subset(select = c("id", "cohort", "state", "txgroup", "time", "month", 
                      "year", "treated", "txtime", "post"))
  
    # Remove excess monthly observations if over monthly_obs (above code
  # automatically creates full years)
  cohort <- subset(cohort, time < start_time + Tobs)
  rownames(cohort) <- 1:nrow(cohort)
  
  return(list("cohort" = cohort,
              "cohort_name" = cohort_name,
              "nperstate" = nperstate,
              "nctrlstates" = nctrlstates,
              "Tobs" = Tobs,
              "Tpost" = Tpost,
              "start_time" = start_time,
              "tStar" = start_time + (Tobs - Tpost)))
}
