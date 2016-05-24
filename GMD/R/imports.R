
## ************************************************************************
## Functions in this file are copied or modified from those in
## other R packages, as stated where appropriate.
## All authorship credit goes to the original author(s).
## 
## Note: A dot is prefixed to the original function names for internal usage.
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Aug 25 14:34:39 EDT 2014 -0400 (Week 34)
## 
## 
## Reference: 
## 
## 
## ************************************************************************



##' A copy of gtools::invalid
##' 
##' see \code{invalid} in package:gtools for details
##' @title Test if a value is missing, empty, or contains only NA or NULL values
##' @param x value to be tested
.invalid <- 
  function(x) 
{
  if (missing(x) || is.null(x) || length(x) == 0) 
    return(TRUE)
  if (is.list(x)) 
    return(all(sapply(x, .invalid)))
  else if (is.vector(x)) 
    return(all(is.na(x)))
  else return(FALSE)
}





##' A copy of wq:::ts2df; see \code{ts2df} in package:wq for details
##' 
##' see \code{ts2df} in package:wq for details.
##' Note: wq_0.3-4 asks for R (>= 2.12.0); but GMD supports R (>= 2.9.0).
##' @title Convert time series to data frame
##' @param x monthly time series vector
##' @param mon1 starting month number, i.e., first column of the data frame
##' @param addYr rows are normally labelled with the year of the starting
##' month, but \code{addYr = TRUE} will add 1 to this year number
##' @param omit if \code{TRUE}, then rows with any \code{NA} will be removed.
ts2df <- 
function(x, mon1 = 1, addYr = FALSE, omit = FALSE) 
{
    if (!is(x, "ts") || is(x, "mts") || !identical(frequency(x), 
        12)) 
        stop("x must be a monthly 'ts' vector")
    if (!mon1 %in% 1:12) 
        stop("mon1 must be between 1 and 12")
    x1 <- window(x, start = c(start(x)[1] - 1, mon1), end = c(end(x)[1] + 
        1, ifelse(mon1 == 1, 12, mon1 - 1)), extend = TRUE)
    d1 <- as.data.frame(matrix(x1, byrow = TRUE, ncol = 12))
    colnames(d1) <- if (mon1 == 1) 
        month.abb
    else month.abb[c(mon1:12, 1:(mon1 - 1))]
    rownames(d1) <- (start(x1)[1] + addYr):(start(x1)[1] + nrow(d1) - 
        1 + addYr)
    d1 = d1[apply(d1, 1, function(x) !all(is.na(x))), ]
    if (omit) 
        d1 = na.omit(d1)
    d1
}



##' Turns a possibly relative file path absolute, performing tilde expansion if necessary.
##'
##' A copy of `tools:::file_path_as_absolute'.
##' @title Turns a Possibly Relative File Path Absolute
##' @return An absolute path
.file_path_as_absolute <- 
  function (x) 
{
  if (length(x) != 1L) 
    stop("'x' must be a single character string")
  if (!file.exists(epath <- path.expand(x))) 
    stop(gettextf("file '%s' does not exist", x), domain = NA)
  cwd <- getwd()
  if (is.null(cwd)) 
    stop("current working directory cannot be ascertained")
  on.exit(setwd(cwd))
  if (file_test("-d", epath)) {
    setwd(epath)
    getwd()
  }
  else {
    setwd(dirname(epath))
    file.path(sub("/$", "", getwd()), basename(epath))
  }
}
