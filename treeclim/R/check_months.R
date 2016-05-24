##' Check if all months are between -12 and 12 (previous and current december)
##' 
##' Make sure that all the months specified in the (probably complex)
##' month specification to dcc are correct months, in the sense that
##' they range between -12 and 12. Further, return the earliest month
##' in the complete selection to decide if we need previous year
##' climate later on.
##' @param selection the month selection as vector, list, or nested
##' list
##' @return a list with the earliest month and the result of the check
##' (logical)
##' @keywords internal
check_months <- function(selection) {
  suppressWarnings(
    try(
      .months <- as.numeric(
        as.vector(
          unlist(
            selection
            )
          )
        )
      )
    )

  out <- list()

  if (length(.months) == 0) {
    
    out$check <- TRUE
    out$minmonth <- -6
    
  } else {
    
    .months <- na.omit(.months)
    if (any(.months < 0)) {
      out$minmonth <- max(.months[which(.months < 0)])
    } else {
      out$minmonth <- min(.months)
    }
    if (!any(.months < -12) & !any(.months > 12)) {
      out$check <- TRUE
    } else {
      out$check <- FALSE
    }
  }
  
  out
}
