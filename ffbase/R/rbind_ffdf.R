rbind.ffdf <- function(..., deparse.level=1){
  a <- list(...)
  x <- ff::clone.ffdf(a[[1]])
  for (l in tail(a, -1)){
    x <- ffdfappend(x, l)
  }
  x
}


#' rbind for ffdf where missing columns are added if not available in one of the ffdf objects
#'
#' rbind for ffdf where missing columns are added if not available in one of the ffdf objects. \cr
#' Similarly as rbind.fill but for ffdf objects
#'
#' @example ../examples/rbindfill.R
#' @param ... 2 or more ffdf objects
#' @param clone logical, indicating to clone the first ffdf object in ... or not before appending
#' the other objects. Defaults to TRUE.
#' @return 
#' an ffdf where the ffdf objects are rbind-ed together. Missing columns in either one of the passed
#' ffdf objects are set to NA values.
#' @export
ffdfrbind.fill <- function(..., clone=TRUE){
  x <- list(...)
  stopifnot(all(sapply(x, FUN=function(x) inherits(x, "ffdf"))))
  columns <- lapply(x, FUN=function(x) colnames(x))
  columns <- do.call(c, columns)
  columns <- unique(columns)
  for(element in 1:length(x)){
    missingcolumns <- setdiff(columns, colnames(x[[element]]))
    for(missingcolumn in missingcolumns){
      x[[element]][[missingcolumn]] <- ff(NA, vmode = "logical", length = nrow(x[[element]]))
    }
  }
  if(clone){
    result <- ff::clone(x[[1]][columns])
  }else{
    result <- x[[1]][columns]
  }
  for (l in tail(x, -1)) {
    result <- ffdfappend(result[columns], l[columns], recode=TRUE)
  }
  result
}

#rbind(as.ffdf(iris), as.ffdf(iris))
#ffdfrbind.fill(as.ffdf(iris), as.ffdf(iris[, 1:3]))
