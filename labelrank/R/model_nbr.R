#' Weights for timing
#' 
#' This function calculates the diminishing weights for label ranking probabilities in case of timing nature of rankings. 
#' @param x a scalar of timing periods.
#' @param n is a parameter of 'memory' of how fast the past gets forgotten.  
#' @return a vector of values.
#' 
#' @details
#' Sometimes rankings have a 'timing' component (for example, weekly sport teams standing) and a recent event can be more important than the past. The model can take advantage of this difference in importance by weighting the ranking probabilities. The weights are calculated using an exponential function . In case of \code{n=1}, weights are a unitary vector of length \code{n}; thus, no time nature in rankings.
time_weights <- function(x, n) {
  n ^ ( ( 1:x) / x - 1)
}

#' A naive Bayes label ranking model
#' 
#' This is an auxiliary function to build a necessary inputs to predict rankings. 
#'  
#' @param x is  \code{n x p} matrix of \code{n} observations and \code{p} training attributes and can have continuous or nominal values.
#' @param y is \code{n x j} matrix of training rankings (permutations). 
#' @param n is a parameter of 'memory'; that is, how fast past gets forgotten. (see details of \link{time_weights}).
#' @return a list of size two: prior and conditional label ranking probabilities.
#' @importFrom stats cor weighted.mean
model_nbr <- function(x, y, n=1) {
    corr <- (cor(t(y[1:nrow(x),]), use = "p") + 1) / 2
    w <- time_weights(nrow(x),n)
    priors <- sapply(1:ncol(corr), function(r) {
        weighted.mean(corr[r, ], w, na.rm = T)
    })
    v <- 1:ncol(x)
    if (all(is.numeric(x))) {
        mu <- w * t(sapply(1:nrow(corr), function(r) {
            sapply(v, function(i) {
                sum(x[, i] * corr[r, ], na.rm = T) / sum(corr[r, ], na.rm = T)
            })
        }))
        sigma <- w * t(sapply(1:nrow(corr), function(r) {
            sapply(v, function(i) {
                sqrt(sum(corr[r, ] * (x[, i] - mu[r, i]) ^ 2, na.rm = T) / sum(corr[r, ], na.rm = T))
            })
        }))
        conditionals <- list(mean = mu, sdev = sigma)
        list(priors = priors, cond = conditionals)
    } else {
        conditionals <- t(do.call(cbind, lapply(v, function(a) {
            sapply(unique(x[, a]), function(value) {
                sapply(1:nrow(x), function(r) {
                  sel.rows <- which(x[, a] == value)
                  sum(corr[sel.rows, r], na.rm = T) / sum(corr[r, ], na.rm = T)
                })
            })
        })))
        list(priors = priors, cond = conditionals)
    }
}
