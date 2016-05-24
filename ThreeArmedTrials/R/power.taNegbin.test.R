#' @title Power related calcuations for three-armed clinical trials with negative binomial distributed endpoints 
#' @description Compute power, sample size, or level of significance for Wald-type test for non-inferiority or superiority of the experimental treatment versus reference treatment with respect to placebo.
#' @details If the individual group sample sizes, i.e. \code{n*allocation} are not natural number, the parameters \emph{n} and \emph{allocation} will be re-calculated.
#' @param rateExp A numeric value specifying the rate of the experimental treatment group in the alternative hypothesis
#' @param rateRef A numeric value specifying the rate of the reference treatment group in the alternative hypothesis
#' @param ratePla A numeric value specifying the rate of the placebo treatment group in the alternative hypothesis
#' @param shape A numeric value specifying the shape parameter
#' @param Delta A numeric value specifying the non-inferiority or superiority margin. Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param sig.level A numeric value specifying the significance level (type I error probability)
#' @param power A numeric value specifying the target power (1 - type II error probability)
#' @param n The total sample size. Needs to be at least 7. 
#' @param allocation A (non-empty) vector specifying the sample size allocation (nExp/n, nRef/n, nPla/n)
#' @param type A character string determing how the variance for the Wald-type test statistic is estimated, must be \emph{restricted}, or \emph{unrestricted} 
#' @return A list with class "power.htest" containing the following components:
#' \item{n}{The total sample size}
#' \item{power}{A numeric value specifying the target power}
#' \item{Delta}{A numeric value specifying the non-inferiority or superiority margin. }
#' \item{sig.level}{A character string specifying the significance level}
#' \item{type}{A character string indicating what type of Wald-type test will be performed}
#' \item{allocation}{A vector with the sample size allocation (nExp/n, nRef/n, nPla/n)}
#' \item{sig.level}{The significance level (Type I error probability)}
#' \item{nExp}{A numeric value specifying the number of sample in the experimental treatment group}
#' \item{nRef}{A numeric value specifying the number of sample in the reference treatment group}
#' \item{nPla}{A numeric value specifying the number of sample in the placebo treatment group}
#' @examples 
#' # Example for type = 'unrestricted': Calculation of n, power, and sig.level. 
#' # Expect n=1038, power=0.8, sig.level=0.025, respectively
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  sig.level = 0.025, power = 0.8, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3))$n
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  sig.level = 0.025, n = 1038, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3))$power
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  power = 0.8007362, n = 1038, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3))$sig.level
#'  
#' # Example for type = 'restricted' calculation of n, power, and sig.level. 
#' # Expect n=1092, power=0.8, sig.level=0.025
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  sig.level = 0.025, power = 0.8, type = 'restricted', allocation = c(1/3, 1/3, 1/3))$n
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  sig.level = 0.025, n = 1092, type = 'restricted', allocation = c(1/3, 1/3, 1/3))$power
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  n = 1092, power = 0.8008113, type = 'restricted', allocation = c(1/3, 1/3, 1/3))$sig.level
#' 
#' # Example for recalculation of 'allocation' and 'n'
#' power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  n = 1001, power = 0.8, allocation = c(0.25, 0.5, 0.25))
#' @export
#' @keywords power NegativeBinomial
power.taNegbin.test <- function(rateExp, rateRef, ratePla, shape, Delta, sig.level = NULL, power = NULL, n = NULL, type = c('restricted', 'unrestricted'), allocation = c(1/3, 1/3, 1/3)){
  
  if( sum(sapply(list(n, power, sig.level), is.null)) !=  1 ) 
    stop("Exactly one of 'n', 'power', and 'sig.level' must be NULL.")
  
  type <- match.arg(type)
  if( missing(rateExp) ){stop("'rateExp' is missing.")}
  if( missing(rateRef) ){stop("'rateRef' is missing.")}
  if( missing(ratePla) ){stop("'ratePla' is missing.")}
  if( missing(shape) ){stop("'shape' is missing.")}
  if( missing(Delta) ){stop("'Delta' is missing.")}
  
  if( !is.numeric(c(rateExp, rateRef, ratePla, shape)) || any(0 >= rateExp, 0 >= rateRef,  0 >= ratePla, 0 >= shape) ){
    stop("'rateExp',  'rateRef', 'ratePla', and 'shape' must not be larger than 0.")
  }
  
  # Calculate effect size and check if parameters are located in the alternative
  effect <- rateExp - Delta * rateRef + (Delta - 1) * ratePla  
  if( effect >= 0 ){
    stop('Parameter vector is not located in the alternative.')
  }

  if( !is.null(sig.level) && !is.numeric(sig.level) || any(0 >= sig.level | sig.level >= 1) ){
    stop("'sig.level' must be numeric in [0, 1]")
  }
  
  if( !is.null(power) && !is.numeric(power) || any(0 >= power | power >= 1) ){
    stop("'power' must be numeric in [0, 1]")
  }
  
  if( !is.null(Delta) && !is.numeric(Delta) || (0 >= Delta) ){
    stop("'Delta' must be larger than 0.")
  }

  if( !is.null(n) && !is.numeric(n) || any(6 >= n) ){
    stop("'n' must be larger than 6.")
  }

  if( !is.null(allocation) && !is.numeric(allocation) || any(allocation>=1,allocation<=0, sum(allocation)!=1, length(allocation)!=3) ){
    stop("'allocation' must not have length 3, sum up to 1, and have only entries between 0 and 1.")
  }
  
  # initialize note
  note <- NULL
  
  # Adjust Sample Size Allocation
  w <- allocation
  if( !is.null(n) ){
    nExp <- ceiling(n * w[1])
    nRef <- ceiling(n * w[2])
    nPla <- ceiling(n * w[3])
    if( any(c(nExp, nRef, nPla) / sum(c(nExp, nRef, nPla)) != w) ){
      n <-  sum(c(nExp, nRef, nPla))
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' and 'n' have been recalculated." 
    }
    else if( sum(c(nExp, nRef, nPla)) != n ){
      n <-  sum(c(nExp, nRef, nPla))
      note <- "'n' have been recalculated." 
    }
  }    
  
  # Calculate variances
  VarUnres <- rateExp * (1+rateExp*shape) / w[1] + Delta^2 * (rateRef * (1+rateRef*shape)) / w[2] + (1-Delta)^2 * (ratePla * (1+ratePla*shape)) / w[3] 
  VarRes <- taNegbin.LimitRestMLE(rateExp1=rateExp, rateRef1=rateRef, ratePla1=ratePla, shape1=shape, Delta=Delta, allocation = allocation)$sigma2.rest

  # Define 'method' for output
  switch(type,
         unrestricted = ( method <- 'Wald-type test with unrestricted variance estimation' ),
         restricted = ( method <- 'Wald-type test with restricted maximum-likelihood variance estimation' )
  )
  
  # Calculate missing parameter
  if( is.null(n) ){
    switch(type,
           unrestricted =  (n <- ceiling( (qnorm(1-sig.level) + qnorm(power))^2 * VarUnres / effect^2) ),
           restricted = (n <- ceiling( (qnorm(1-sig.level)*sqrt(VarRes)/ sqrt(VarUnres) + qnorm(power))^2 * VarUnres / effect^2) )
    )
    nExp <- round(n * w[1])
    nRef <- round(n * w[2])
    nPla <- round(n * w[3])
    n <- nExp + nRef + nPla
    if( any(c(nExp, nRef, nPla) / n != w) ){
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' has been recalculated." 
    }
  }  
  if( is.null(power) ){
    switch(type,
           unrestricted =  (power <- pnorm(qnorm(sig.level) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1) ),
           restricted = (power <- pnorm(qnorm(sig.level)*sqrt(VarRes)/ sqrt(VarUnres) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1))
    )
  }
  if( is.null(sig.level) ){
    switch(type,
           unrestricted =  (sig.level <- pnorm(qnorm(power) + sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1) ),
           restricted = (sig.level <- pnorm(qnorm(power)*sqrt(VarUnres)/ sqrt(VarRes) + sqrt(n) * effect / sqrt(VarRes), mean = 0, sd = 1))
    )
  }
  
  structure(list(method = method,
                 rateExp = rateExp,
                 rateRef = rateRef,
                 ratePla = ratePla,
                 shape = shape,
                 n = n,
                 sig.level = sig.level,
                 power = power,
                 type = type,
                 Delta = Delta,
                 allocation = w,
                 nExp = nExp, 
                 nRef = nRef,
                 nPla = nPla,
                 note = note),
            class = "power.htest")
}


