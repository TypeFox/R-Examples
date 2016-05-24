###############################################
# --------------------------------------------#
# Functions for the doubly t-distribution     #
# --------------------------------------------#
###############################################

#' Approximation of the distribution function of the doubly noncentral t-distribution 
#'
#' Computes the value of the distribution function of the doubly noncentral t-distribution at \code{x}.
#'
#' @inheritParams overview
#' 
#' @return Distribution value of the doubly noncentral t-distribution at \code{x}.
#'
#' @export
doublyT <- function(x, df, delta, lambda, lb = 0, ub) {
  values <- lb:ub
  suppressWarnings(sum(dpois(values, lambda/2)*pt(x*(1+2*values/df)^(0.5), df+2*values, delta)))
}

#' Calculation of the NCPs of each randomization sequence for the doubly noncentral t-distribution
#'
#' Computes the noncentraility parameters delta and lambda for the doubly noncentral t-distribution of each randomization sequence.
#'
#' @param randSeq object of the class randSeq.
#' @param bias object of the class bias.
#' @param endp object of the class endpoint.
#' 
#' @examples 
#' myPar <- crPar(4)
#' M <- getAllSeq(myPar)
#' cs <- selBias("CS", 1, "exact")
#' endp <- normEndp(mu = c(0, 0), sigma = c(1, 1))
#' genNcps(M, cs, endp)
#' 
#' @return matrix containing the noncentrality parameters delta and lambda of all randomization sequences.
#'
#' @export
genNcps <- function(randSeq, bias, endp) {
  stopifnot(is(randSeq, "randSeq"), randSeq@K == 2,
            is(bias, "issue"), is(endp, "endpoint"), sum(duplicated(endp@sigma)) == 1)
  # all randomization sequences
  allRandSeq <- randSeq@M
  # all expectation values in form of a matrix
  allExp <- getExpectation(randSeq, bias, endp)
  # number of randomization sequences
  r <- nrow(allRandSeq)

  P <- lapply(1:r, function(i) {
      randSeq <- allRandSeq[i, ]
      exp <- allExp[i, ]
      splitGroups <- split(exp, randSeq)
      # average expectations in both groups
      avExp <- lapply(splitGroups, mean)
      # variance of the expectations in both groups
      varExp <- lapply(splitGroups, function(x) sum((x-mean(x))^2))
      # number of assigend patients to both of the groups
      numAssPat <- lapply(splitGroups, length)
      # defining P
      P <- matrix(0, nrow = 1, ncol = 5)
      # updating the matrix P
      P[1, 3] <- do.call("sum", numAssPat)
      P[1, 4] <- ifelse(!is.null(numAssPat$"0"), numAssPat$"0", 0)
      P[1, 5] <- ifelse(!is.null(numAssPat$"1"), numAssPat$"1", 0)
      
      if (!is.null(numAssPat$"0") && !is.null(numAssPat$"1")) {
        P[1, 1] <- 1/endp@sigma[1] * sqrt((P[1, 4]*P[1, 5]) / (P[1, 4]+P[1, 5])) * (avExp$"0" - avExp$"1")
        P[1, 2] <- 1/endp@sigma[1]^2 * (varExp$"0" + varExp$"1")
        return(P)
      } else {
        P[1, 1] <- 0
        P[1, 2] <- sum((exp - mean(exp))^2)/endp@sigma[1]^2
        return(P)
      } 
    }
  )
  P <- do.call("rbind", P)
  colnames(P) <- c("delta", "lambda", "N", "nA", "nB")
  P
}


#' Calculation of the biased type-one-error (resp. power) of Student`s t-test
#'
#' Computes the biased type-one-error (resp. power) of Student`ts t-test due to shifts in the expectation vectors 
#' in both treatment groups.
#'
#' @param randSeq object of the class randSeq.
#' @param bias object of the class bias.
#' @param endp object of the class endpoint.
#' 
#' @examples 
#' myPar <- crPar(4)
#' M <- getAllSeq(myPar)
#' cs <- selBias("CS", 1, "exact")
#' endp <- normEndp(mu = c(0, 0), sigma = c(1, 1))
#' doublyTValues(M, cs, endp)
#' 
#' @return the biased type-one-error (resp. power) of all randomization sequences.
#'
#' @export
doublyTValues <- function(randSeq, bias, endp) {
  stopifnot(is(randSeq, "randSeq"), randSeq@K == 2,
            is(bias, "issue"), is(endp, "endpoint"), sum(duplicated(endp@sigma)) == 1)   
  # calculation of of the matrix containing the ncp
  ncps <- genNcps(randSeq, bias, endp)
  # variable for the defined alpha quantile
  alpha <- bias@alpha
  # function for calculating the p values of the singular randomization sequences
  p.value <- sapply(1:nrow(ncps), function(i) {
      x <- ncps[i, ]
      # return zero if in one treatment group was no observation
      if( x[4] == 0 || x[5] == 0)
        return(0)
      
      # lower boundary
      lb <- max(floor(x[2]/2 - qpois(.995, x[2]/2)), 0)
      # upper boundary
      ub <- as.vector(ceiling(x[2]/2 + qpois(.995, x[2]/2)))
      # degrees of freedom
      df <- as.vector(x[3]) - 2
      # t quantiles
      tQuantLow <- qt(alpha/2, df)
      tQuantUpper <- -tQuantLow 
      p.value.less <- doublyT(tQuantLow , df, as.vector(x[1]), as.vector(x[2]), lb, ub)
      p.value.greater <- 1 - doublyT(tQuantUpper, df, as.vector(x[1]), as.vector(x[2]), lb, ub)
      p.value <- p.value.less + p.value.greater
      return(p.value)
    }
  )
  p.value
}

