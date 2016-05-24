#' First order features
#'
#' @param data Numeric 2D matridata.
#' @references \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102107#s5} 
#' @name first_order_features
NULL
#> NULL

#' @describeIn first_order_features Energy (ASM)
#' 
calc_energy <- function(data){
  #TODO: Add dim check for 2D vs 3D
  return(sum(as.numeric(data)*as.numeric(data), na.rm=TRUE))
}

#' @describeIn first_order_features Entropy
#' @param base The base for which the logarithm is calculate
#' @param nbins The number of bins the histogram is discretized into
calc_entropy <- function(data, base=2, nbins=length(unique(c(data)))){
  # Break data into a hist
    im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data[!is.na(data)])
  
  #Logs cannot take 0 values, so let = 0 if no value
  entropy_vals <- vapply(cuts, function(data) ifelse(data != 0, data*logb(data, base=base), 0), FUN.VALUE = 1)
  return(-1*sum(entropy_vals))
}

#' @describeIn first_order_features Kurtosis
#' 
calc_kurtosis <- function(data){
  n <- length(data[!is.na(data)])
  data <- data - mean(data, na.rm=TRUE)
  r <- n * sum(data^4, na.rm=TRUE) / (sum(data^2, na.rm=TRUE)^2)
  return(r * (1 - 1/n)^2 - 3)
}

#' @describeIn first_order_features Mean Deviation
#' 
calc_meanDeviation <- function(data){
  scale <- 1/prod(dim(data))
  mu <- mean(data, na.rm=TRUE)
  return(scale * sum(abs(data - mu), na.rm=TRUE))
}

#' @describeIn first_order_features Skewness
#' 
calc_skewness <- function (data){
  
  data <- data[!is.na(data)]

  return(sum((data - mean(data))^3)/(length(data) * sd(data)^3))
}

#' @describeIn first_order_features Uniformity
#'
calc_uniformity <- function(data, nbins=length(unique(c(data)))){
  # Break data into a hist
  data <- data[!is.na(data)]
  im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data)
  function_vals <- vapply(cuts, function(data) data^2, FUN.VALUE = 1)
  return(sum(function_vals))
}

#' @describeIn first_order_features Mean
#'
calc_mean <- function(data) mean(data, na.rm=TRUE)

#' @describeIn first_order_features Median
#'
calc_median <- function(data) median(data, na.rm=TRUE)

#' @describeIn first_order_features Maximum Value
#'
calc_max <- function(data) max(data, na.rm=TRUE)

#' @describeIn first_order_features Minimum Value
#'
calc_min <- function(data) min(data, na.rm=TRUE)

#' @describeIn first_order_features Variance
#'
calc_variance <- function(data) var(c(data), na.rm=TRUE)

#' @describeIn first_order_features Root Mean Squared
#'
calc_RMS <- function(data) sqrt(mean(data^2, na.rm=TRUE))

#' @describeIn first_order_features Standard Deviation
#'
calc_sd <- function(data) sd(data, na.rm=TRUE)
