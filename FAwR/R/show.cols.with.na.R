show.cols.with.na <-
function(x) {
  if (class(x) != "data.frame") 
    stop("x must be a dataframe.\n")
  missing.by.column <- apply(is.na(x), 2, sum)
  if (sum(missing.by.column) == 0) {
    cat("No missing values.\n")
  } else {
    missing <- which(missing.by.column > 0)
    return(missing.by.column[missing])
  }
}

