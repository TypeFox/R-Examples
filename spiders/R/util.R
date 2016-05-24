##' sum specified columns by time
##'
##' @param data a dataframe
##' @param vars column variables in data to sum over
getTimeCounts <- function(data, vars) {
    ## if ( missing(by) ) {
    ##     by <- 'time'
    ## } else {
    ##     by <- c('time', by)
    ## }
    ddply(data, "time", colwise(sum, vars))
}

##' checks convergence of the parameters for the estimation functions
##'
##' @param theta an object, convertible to a matrix, of current parameter values
##' @param theta_old an object, convertible to a matrix, of old parameter values
##' @param eps tolerance to determine convergence
converged <- function(theta, theta_old, eps=1e-5) {
    isTRUE(all.equal(theta, theta_old, tolerance=eps, check.names=F, check.attr=F))
}

## sum over species to get a vector of values for each time period
sumSp <- function(mat) {
    matrix(rowSums(mat), nrow=nrow(mat))
}

## sum over times to get a vector of values for each species
sumT <- function(mat) {
    matrix(colSums(mat), ncol=ncol(mat))
}

sumST <- function(mat) {
    sum(mat)
}

## colors stolen from http://geography.uoregon.edu/datagraphics/color_scales.htm
cols <- c(orange1="#FFBF80", orange2="#FF8000",
          yellow1 = "#FFFF99", yellow2="#FFFF33",
          green1 = "#B2FF8C", green2="#33FF00",
          blue1 = "#A6EDFF", blue2="#1AB2FF",
          purple1 = "#CCBFFF", puple2="#664CFF",
          red1 = "#FF99BF", red2="#E61A33")
