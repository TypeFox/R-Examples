# Wrapper functions for applying linear interpolation, and correlation and 
# kolmogorov-smirnov tests to multiple variables in separate datasets.

#' Interpolate multiple columns of a data.frame
#' 
#' @param data1 a \code{data.frame} with columns for interpolation
#' @param n the number of points to interpolate
#' @param method the method to interpolate by, \code{"linear"} or 
#'   \code{"spline"}
#' @return Returns a \code{data.frame} of all \code{data1} variables
#'   interpolated to \code{n} rows.
#' @seealso \code{\link{approx}} is the underlying function and has more
#'   options
#' @importFrom stats approx spline
#' @export
#' 
approxm <- function(data1, n, method = "linear") {
  newdata <- data.frame(matrix(nrow = n, 
                               ncol = ncol(data1)))
  
  colnames(newdata) <- names(data1)
  
  for(i in names(data1)) {
    newdata[, i] <- if(method == "linear") {
      approx(data1[, i], n = n)$y
    } else if(method == "spline") {
      spline(data1[, i], n = n)$y
    }
  }
  
  return(newdata)
}


#' Kolmogorov-Smirnov test for multiple variables - of the same name - in 
#' separate data.frames
#' 
#' @param data1 a \code{data.frame} with N variables
#' @param data2 a \code{data.frame} with the same N variables as data1
#' @param vars a \code{vector} of the N variable names
#'   
#' @return Returns a \code{data.frame} with Kolmogorov-Smivnov test data for N
#'   variables
#' @seealso \code{\link{ks.test}} is the underlying function and has more
#'   options.
#' @importFrom stats ks.test
#'   
#' @export
#' 
ks.testm <- function(data1, data2, vars) {
  a <- data.frame(matrix(nrow = length(vars), 
                         ncol = 2))
  
  for(i in 1:length(vars)) {
    a[i, 1:2] <- unlist(ks.test(data1[, vars[i]],
                                data2[, vars[i]]))[1:2]
  }
  
  a <- cbind(vars, a)
  names(a) <- c("vars", "D", "P")
  
  return(a)
}

#' Correlation test for multiple variables - of the same name - in separate 
#' data.frames
#' 
#' @param data1 a \code{data.frame} with N variables
#' @param data2 a \code{data.frame} with the same N variables as \code{data1}
#' @param method a correlation method, either \code{method = "pearson"} or 
#'   \code{"spearman"}
#'   
#' @return Returns a \code{data.frame} of correlation test data for N variables
#' @seealso \code{\link{cor.test}} is the underlying function and has more
#'   options
#' @importFrom stats cor.test
#' @export
#' 
cor.testm <- function(data1, data2, method) {
  newdata <- list()
  x = 0
  for(i in names(data1)) {
    x = x + 1
    newdata[[x]] <- cor.test(data1[, i], 
                             data2[, i], 
                             method = method)
  }
  
  # Reshape data
  newdata <- data.frame(t(sapply(newdata, c)))
  newdata <- apply(newdata, 2, unlist)
  
  # Split conf.int (double length vector) and remove original
  if(method == "pearson") {
    ci     <- split(newdata$conf.int, 1:2)
    ci.min <- ci[1]
    ci.max <- ci[2]
    newdata$conf.int <- NULL
  }
  
  # Add and rearrange variables for readability
  newdata <- data.frame(var    = names(data1),
                        cor    = newdata[["estimate"]],
                        ci.min = ifelse(method == "pearson", ci.min, NA),
                        ci.max = ifelse(method == "pearson", ci.max, NA),
                        p      = newdata[["p.value"]],
                        stat   = newdata[["statistic"]],
                        df     = ifelse(method == "pearson", 
                                        newdata[["parameter"]], NA),
                        alt    = newdata[["alternative"]],
                        method = method)
  names(newdata)[3:4] <- c("ci.min", "ci.max")
  return(newdata)
}
