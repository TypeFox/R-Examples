##' Get positions and names for months
##' 
##' Get positions for supplied numeric month ids in the 2-year span
##' under consideration, and connect to abbreviated month names for
##' previous and current year.
##' 
##' @param month range or set of numeric month ids
##' @return a list of positions and names
##' @keywords manip internal
format_month <- function(month) {
  ## 2 possibilities: continuous or non-continuous. if it is not
  ## continuous, it will either be correctly specified by the user, or
  ## it will be supplied correctly by exclude_from. it is not possible
  ## to mix ranges through zero with other specifications (e.g.,
  ## c(-4:3, 1, 5) is not possible).
  if (is_continuous(month)) {
    month <- correct_continuous(month)
  } else {
    if (any(month == 0)) {
      stop("It is not possible to mix ranges through zero with other specifications.")
    }
  }
  
  ## month ids and names
  month_ids <- c(-1:-12, 1:12)
  ## check if months are between -1 and 12
  if (any(!is.element(month, month_ids))) {
    stop("Month specification must be within previous january (-1) and current december (12).")
  }

  mmonth <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
              "sep", "oct", "nov", "dec")
  ucmonth <- toupper(mmonth)
  
  two_years <- paste(
    c(
      rep("prev.", 12),
      rep("curr.", 12)),
    rep(mmonth,
        2), sep = "")
  single <- c(mmonth, ucmonth)
  .month <- list()
  .month$match <- match(month, month_ids)
  .month$names <- two_years[.month$match]
  .month$single <- single[.month$match]
  .month
}
