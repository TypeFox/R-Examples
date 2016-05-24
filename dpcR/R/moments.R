#' Calculate Moments of Poisson Distribution
#' 
#' Computes moments of a Poisson distribution. The calculations are based on 
#' values of positive and total partitions or the theoretical lambda value.
#' 
#' 
#' @name moments-methods
#' @aliases moments-methods moments,dpcr-method moments,numeric-method 
#' moments,matrix-method
#' @docType methods
#' @param input a vector with two elements (the first element is treated as 
#' a number of positive partitions and the second as a number of total partitions) 
#' or a matrix with two columns (first columns contains numbers of 
#' positive partitions and the second total numbers of total partitions)
#' or an object of class \code{\linkS4class{dpcr}} .
#' @return A data frame with four columns: name of the experiment, 
#' name of the replicate, method of computation (theoretical or empirical), 
#' name of the moment and the value of the moment. The theoretical moments are 
#' computed using the lambda value and the empirical using the sample values.
#' @note Four first moments of a Poisson distribution.
#' 
#' Mean : \eqn{\lambda}{lambda}.
#' 
#' Variance: \eqn{\lambda}{lambda}.
#' 
#' Skewness: \eqn{\sqrt{\lambda}}{lambda^(-0.5)}.
#' 
#' Kurtosis: \eqn{\frac{1}{\lambda}}{lambda^(-1)}.
#' @author Michal Burdukiewicz.
#' @keywords moments mean variance skewness kurtosis
#' @export
#' @examples
#' # moments for 100 positive partitions of 765 total partitions
#' moments(c(100, 765))
#' 
#' # calculate moments for an array digital PCR
#' moments(six_panels)
#' 
moments <- function (input) {
  stop("Wrong class of 'input'")
}

setMethod("moments", signature(input = "numeric"), function(input) {
  if (length(input) != 2) 
    stop("Input must contain two values: number of positive paritions and total number 
         of partitions.") 
  

  calc_moms(data.frame(experiment = "Unknown", replicate = "Unknown",
                       k = input[1], n = input[2]))
})

setMethod("moments", signature(input = "matrix"), function(input) {
  if (ncol(input) != 2) 
    stop("Input matrix must contain two columns: a number of positive paritions and 
         a total number of partitions.") 

  calc_moms(data.frame(experiment = 1L:nrow(input), replicate = "Unknown",
                       k = input[, 1], n = input[, 2]))
})


setMethod("moments", signature(input = "dpcr"), function(input) {
  summ <- summary(input, print = FALSE)[["summary"]]
  calc_moms(summ[summ[["method"]] == "bhat", c("experiment", "replicate", "k", "n")])
})

# kn is a data.frame of length 4, with specified experiment name, replicate id,
# number of positive partitions and total number of partitions
calc_moms <- function(kn) {
  do.call(rbind, lapply(1L:nrow(kn), function(single_row) 
    cbind(experiment = rep(kn[single_row, "experiment"], 8),
          replicate = rep(kn[single_row, "replicate"], 8),
          rbind(data.frame(method = rep("theoretical", 4), 
                           moment = c("mean", "variance", "skewness", "kurtosis"),
                           value = moms(kn[single_row, "k"], kn[single_row, "n"])),
                data.frame(method = rep("empirical", 4), 
                           moment = c("mean", "variance", "skewness", "kurtosis"),
                           value = empir_moms(kn[single_row, "k"], kn[single_row, "n"]))
          )
    )
  ))
}


#first four moments of distribution
#kn_df - data frame, row is run, first column is k, second is n
moms <- function(k, n) {
  lambda <- fl(k/n)
  c(lambda, lambda, lambda^(-0.5), 1/lambda)
}

#sample moments
empir_moms <- function(k, n) {
  exp_kn <- c(rep(1, k), rep(0, n - k))
  #res <- c(k/n, k(n - k)*n^-2, skewness(input), kurtosis(input))
  c(mean(exp_kn), var(exp_kn), skewness(exp_kn), kurtosis(exp_kn))
}