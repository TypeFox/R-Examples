# stacks.r

#' Multinomial group testing estimation of stacked genes
#'
#' Assuming qualitative tests are performed on n pools of m seeds, use
#' multinomial group testing to estimate the percent of seeds with single
#' genetic traits and the percentage of seeds with stacked genetic traits.
#'
#' The 'stack2Excel' and 'stack3Excel' functions are simple wrappers that are
#' intended to be called from Excel and should not issue any warnings.
#'
#' @rdname stack2
#' @aliases stack2 stack3 stack2Excel stack3Excel
#' @param n the number of pools
#' @param m the number of seeds in each pool
#' @param nA the number of positive pools for event A only
#' @param nB the number of positive pools for event B only
#' @param nAB the number of positive pools for both A and B
#' @param nC the number of positive pools for event C only
#' @param nAC the number of positive pools for both A and C
#' @param nBC the number of positive pools for both B and C
#' @param nABC the number of positive pools for both A and B and C
#' @param existAB do seeds with a stacked event 'AB' exist?
#' @param existAC do seeds with a stacked event 'AC' exist?
#' @param existBC do seeds with a stacked event 'BC' exist?
#' @param existABC do seeds with a stacked event 'ABC' exist?
#' @param fpr false positive rate (proportion) for detecting GM events
#' @param fnr false negative rate (proportion) for detecting GM events
#' @param check Should simple checks be performed?  Defaults to TRUE
#' @return A data frame with the estimated proportion of seeds for each event,
#' the observed and expected number of positive pools, and whether or not each
#' event can exist.
#' @author Kevin Wright, Jean-Louis Laffont
#' @examples
#' stack2(10, 300, 0, 1, 2)
#' stack3(20,150, 2,2,2,2,2,2,3, existAB="no", fnr=.02, fpr=.02)
#' @export
stack3 <- function(n, m, nA, nB, nC, nAB, nAC, nBC, nABC,
                       existAB="Yes", existAC="Yes", existBC="Yes", existABC="Yes",
                       fpr=0, fnr=0, check=TRUE) {
	# ML estimation of the proportions of GM seeds when THREE events are tested
	# n: Total number of pools
	# m: number of seeds per pool
	# nA, nB, nC: number of positive pools for A, B, C only
	# nAB: number of positive pools for both A and B
  # existAB: does stack AB exist?
  # fnr:   false negative rate (%)
  # fpr:   false positive rate (%)

  if (check && sum(c(nA, nB, nC, nAB, nAC, nBC, nABC)) > n)
    stop("The total number of positive pools must be <=", n)

  # Define indicator vector for which 'theta' to optimize
  ind <- c(FALSE, TRUE, TRUE, TRUE, # 0, A, B, C
           existAB=="Yes", existAC=="Yes", existBC=="Yes", existABC=="Yes")

  N <- sum(c(nA, nB, nC, nAB, nAC, nBC, nABC))
  n0 <- n-N
  nPos <- c(n0, nA, nB, nC, nAB, nAC, nBC, nABC)
	# Initial values and bounds of thetas to be optimized
  init <- (1 - (1 - N/n)^(1/m)) * nPos[ind]/N
  names(init) <- c('0', 'A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC')[ind]
  upper <- rep(1, length(init))
  lower <- rep(1e-7, length(init))  # needs to be larger than ndeps

  cpr = 1-fpr # Correct Positive Rate
  cnr = 1-fnr # Correct Negative Rate
  
  K <- matrix(c(cpr^3,     fnr*cpr^2,   fnr*cpr^2,   fnr*cpr^2,   fnr^2*cpr,   fnr^2*cpr,   fnr^2*cpr,   fnr^3,
                fpr*cpr^2, cnr*cpr^2,   fnr*fpr*cpr, fnr*fpr*cpr, fnr*cnr*cpr, fnr*cnr*cpr, fpr*fnr^2,   cnr*fnr^2,
                fpr*cpr^2, fnr*fpr*cpr, cnr*cpr^2,   fnr*fpr*cpr, fnr*cnr*cpr, fpr*fnr^2,   fnr*cnr*cpr, cnr*fnr^2,
                fpr*cpr^2, fnr*fpr*cpr, fnr*fpr*cpr, cnr*cpr^2,   fpr*fnr^2,   fnr*cnr*cpr, fnr*cnr*cpr, cnr*fnr^2,
                fpr^2*cpr, cnr*fpr*cpr, cnr*fpr*cpr, fnr*fpr^2,   cnr^2*cpr,   fnr*cnr*fpr, fnr*cnr*fpr, cnr^2*fnr,
                fpr^2*cpr, cnr*fpr*cpr, fnr*fpr^2,   cnr*fpr*cpr, fnr*cnr*fpr, cnr^2*cpr,   fnr*cnr*fpr, cnr^2*fnr,
                fpr^2*cpr, fnr*fpr^2,   cnr*fpr*cpr, cnr*fpr*cpr, fnr*cnr*fpr, fnr*cnr*fpr, cnr^2*cpr,   cnr^2*fnr,
                fpr^3,     cnr*fpr^2,   cnr*fpr^2,   cnr*fpr^2,   cnr^2*fpr,   cnr^2*fpr,   cnr^2*fpr,   cnr^3),
              ncol=8, byrow=TRUE)

  oo <- nlminb(init, nll3,
               nPos=nPos, m=m, K=K,
               lower=lower, upper=upper, control = list(rel.tol = 1e-7))
  ok <- oo$convergence==0 && is.finite(oo$objective)

  # Assemble results
  dat <- data.frame(event = c('A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC'),
                    prop = NA,
                    obs = nPos[-1],
                    expect=NA,
                    exists = ind[-1])
  if(ok){
    dat$prop <- rep(0, nrow(dat))
    dat$prop[match(names(oo$par),dat$event)] <- oo$par

    # Expected number of pools
    tA <- dat$prop[1] ; tB <- dat$prop[2] ; tC <- dat$prop[3]
    tAB <- dat$prop[4] ; tAC <- dat$prop[5] ;   tBC <- dat$prop[6]
    tABC <- dat$prop[7]
    p0 <- valid((1 - (tA + tB + tC + tAB + tAC + tBC + tABC))^m)
    pA <- valid((1 - (tB + tC + tAB + tAC + tBC + tABC))^m - p0)
    pB <- valid((1 - (tA + tC + tAB + tAC + tBC + tABC))^m - p0)
    pC <- valid((1 - (tA + tB + tAB + tAC + tBC + tABC))^m - p0)
    pAB <- valid((1 - (tC + tAC + tBC + tABC))^m - p0 - pA - pB)
    pAC <- valid((1 - (tB + tAB + tBC + tABC))^m - p0 - pA - pC)
    pBC <- valid((1 - (tA + tAB + tAC + tABC))^m - p0 - pB - pC)
    pABC <- valid(1 - p0 - pA - pB - pC - pAB - pAC - pBC)
    dat$expect <- n * c(pA, pB, pC, pAB, pAC, pBC, pABC)
  }
  class(dat) <- c("seedstack", class(dat))
  return(dat)
}

# ----------------------------------------------------------------------------

#' #' Ensure probabilities are valid
#'
#' Force calculated probabilities into the range [0,1].
#'
#' Due to floating-point arithmetic, a number that should represent a
#' probability can be calculated as being less than zero or greater than one.
#' This function returns a value that is a valid probability.
#'
#' @param x probability
valid <- function(x){
  # Restrict to 0,1.  Note: pmin/pmax was slow, so was ifelse.
  return(max(0, min(1, x)))
}

# ----------------------------------------------------------------------------


#' @param ... Other arguments passed
#' @rdname stack2
#' @export
stack2Excel <- function(...){
  stack2(..., check=FALSE)
}

#' @rdname stack2
#' @export
stack3Excel <- function(...){
  stack3(..., check=FALSE)
}

# ----------------------------------------------------------------------------

#' @rdname stack2
#' @export
stack2 <- function(n, m, nA, nB, nAB, existAB="Yes", fpr=0, fnr=0,
                   check=TRUE) {
	# Estimation of the proportions of GM seeds when TWO events are tested
	# n: number of pools
	# m: number of seeds per pool
	# nA, nB: number of positive pools for A, B only
	# nAB: number of positive pools for both A and B
  # existAB: does stack AB exist?
  # fnr: false negative rate
  # fpr: false positive rate

  if(check && sum(nA, nB, nAB) > n)
    stop("The total number of positive pools must be <=",n)

  # Indicator vector for which values of theta to optimize
  ind <- c(FALSE, TRUE, TRUE, # 0, A, B
           existAB=="Yes")

  N <- sum(c(nA, nB, nAB))
  n0 <- n-N
  nPos <- c(n0, nA, nB, nAB)
	# Initial values and bounds of thetas to be optimized
  init <- (1 - (1 - N/n)^(1/m)) * nPos[ind]/N
  names(init) <- c('0', 'A', 'B', 'AB')[ind]
  upper <- rep(1, length(init))
  lower <- rep(1e-7, length(init))  # needs to be larger than ndeps

  cpr = 1-fpr # correct positive rate
  cnr = 1-fnr
  K <- matrix(c(cpr^2,   fnr*cpr, fnr*cpr, fnr^2,
                fpr*cpr, cpr*cnr, fpr*fnr, fnr*cnr,
                fpr*cpr, fpr*fnr, cpr*cnr, fnr*cnr,
                fpr^2,   fpr*cnr, fpr*cnr, cnr^2), ncol=4, byrow=TRUE)

  oo <- nlminb(init, nll2,
               nPos=nPos, m=m, K=K,
               lower=lower, upper=upper, control=list(rel.tol=1e-7))
  ok <- oo$convergence==0 && is.finite(oo$objective)

  # Assemble results
  dat <- data.frame(event = c('A', 'B', 'AB'),
                    prop = NA,
                    obs = nPos[-1],
                    expect = NA,
                    exists = ind[-1])
  if(ok){
    dat$prop <- rep(0, nrow(dat))
    dat$prop[match(names(oo$par),dat$event)] <- oo$par

    # Expected number of pools
    tA <- dat$prop[1] ; tB <- dat$prop[2] ; tAB <- dat$prop[3]
    p0 <- (1 - (tA + tB + tAB))^m
    pA <- (1 - (tB + tAB))^m - p0
    pB <- (1 - (tA + tAB))^m - p0
    pAB <- 1 - p0 - pA - pB
    dat$expect <- n * c(pA, pB, pAB)
  }
  class(dat) <- c("seedstack", class(dat))
  return(dat)
}

nll2 <- function(theta, nPos, m, K){
  # Negative loglikelihood
  # theta   = [thetaA, thetaB, thetaAB] (or some subset)
  # nPos    = [n0, nA, nB, nAB]
  # m       = seeds per pool
  # K       = coefficients to fnr/fpr adjustment

  # Optimize with respect to certain theta values
  # Missing names(theta) return NA, which are changed to 0.
  tA <- theta["A"]
  tB <- theta["B"]
  tAB <- theta["AB"]   ; tAB <- ifelse(is.na(tAB), 0, tAB)

  # Switch from theta to p
  p0 <- valid((1 - (tA + tB + tAB))^m)
  pA <- valid((1 - (tB + tAB))^m - p0)
  pB <- valid((1 - (tA + tAB))^m - p0)
  pAB <- valid(1 - p0 - pA - pB)

  # Adjust for fpr/fnr.
  probs <- K %*% c(p0, pA, pB, pAB)
  probs <- pmax(0, probs)
  probs <- pmin(1, probs)

  # ll: n0*log(p0) + nA*log(pA) + nB*log(pB) + nAB*log(pAB)
  if(any(is.nan(probs))) {
    ll <- NA
  } else if (all(probs > 0)) {
    ll <- sum(nPos * log(probs))
  } else ll <- -Inf

  return(-ll)
}

nll3 <- function(theta, nPos, m, K){
  # theta   [thetaA, thetaB, ... thetaABC] (or some subset)
  # nPos    [n0, nA, nB, ... , nABC]
  # m       seeds per pool
  # K       coefficients to fnr/fpr adjustment

  # Optimize with respect to certain theta values
  # Missing names(theta) return NA, which are changed to 0.
  tA <- theta["A"]
  tB <- theta["B"]
  tC <- theta["C"]
  tAB <- theta["AB"]   ; tAB <- ifelse(is.na(tAB), 0, tAB)
  tAC <- theta["AC"]   ; tAC <- ifelse(is.na(tAC), 0, tAC)
  tBC <- theta["BC"]   ; tBC <- ifelse(is.na(tBC), 0, tBC)
  tABC <- theta["ABC"] ; tABC <- ifelse(is.na(tABC), 0, tABC)

  # Switch from theta to p
  p0 <- valid((1 - (tA + tB + tC + tAB + tAC + tBC + tABC))^m)
  pA <- valid((1 - (tB + tC + tAB + tAC + tBC + tABC))^m - p0)
  pB <- valid((1 - (tA + tC + tAB + tAC + tBC + tABC))^m - p0)
  pC <- valid((1 - (tA + tB + tAB + tAC + tBC + tABC))^m - p0)
  pAB <- valid((1 - (tC + tAC + tBC + tABC))^m - p0 - pA - pB)
  pAC <- valid((1 - (tB + tAB + tBC + tABC))^m - p0 - pA - pC)
  pBC <- valid((1 - (tA + tAB + tAC + tABC))^m - p0 - pB - pC)
  pABC <- valid(1 - p0 - pA - pB - pC - pAB - pAC - pBC)

  # Adjust for fpr/fnr.
  probs <- K %*% c(p0, pA, pB, pC, pAB, pAC, pBC, pABC)
  probs <- pmax(0, probs)
  probs <- pmin(1, probs)

  # ll: n0*log(p0) + nA*log(pA) + ... + nABC*log(pABC)
  if(any(is.nan(probs))) {
    ll <- NA
  } else if(all(probs > 0)) {
    ll <- sum(nPos * log(probs))
  } else ll <- -Inf

  return(-ll)
}


# ----------------------------------------------------------------------------


#' Print stack object.
#'
#' Print method for seedstack object.
#'
#' @param x A data frame to print pretty.
#' @rdname stack2
#' @method print seedstack
#' @export
print.seedstack <- function(x, ...) {
  # formatC chokes on c(NA,NA), so use as.numeric first
  x <- data.frame(event=c(as.character(x$event),"Total"),
                    percent=paste(formatC(as.numeric(100*c(x$prop, sum(x$prop))), digits=2, format='f'), "%"),
                    obs=c(prettyNum(x$obs),""),
                    expect=c(formatC(as.numeric(x$expect), digits=1, format='f'),""),
                    exists=c(ifelse(x$exists,"Yes","No"),"")
                    )
  colnames(x) <- c("Event", "Estimate", "Observed", "Expected", "Exists")

  cat("Percent GM seeds for each event:\n")
  print.data.frame(x, row.names=FALSE)
  return()
}

if(FALSE){

stack3(20,150, 2,2,2,2,2,2,3, existAB="no", fnr=.02, fpr=.02) # paper

stack2(20, 150, 8, 5, 1)

stack2(20, 150, 8, 5, 8) # fails in R
stack2(20, 150, 8, 5, 7, check=FALSE) # fails
stack2Excel(20, 150, 8, 5, 7) # fails

stack2(20, m=150, 3, 0, 0) # estimates are all NA
stack2(20, m=150, 10, 0, 0, existAB='no') # two digits in 'observed' column

stack3(20, 150, 1,2,3,1,2,1,1)
stack2(20, 100, 3, 3, 3, existAB='no')
stack2(20, 150, 3, 3, 3, existAB='Yes')

stack2(10, 300, 0, 1, 2) # JLL example 1
stack3(20,150, 2,2,2,2,2,2,3, existAB="no", fnr=.02, fpr=.02) # JLL example 2

stack3(20,150, 2,2,2,2,2,2,3, existAB='no')
stack3(20,150, 2,2,2,2,2,2,3, existAC='no')
stack3(20,150, 2,2,2,2,2,2,3, existBC='no')
stack3(20,150, 2,2,2,2,2,2,3, existABC='no')

stack3(20,150, 2,2,2,2,2,2,9, existABC='no') # should fail in R
stack3(20,150, 2,2,2,2,2,2,9, existABC='no', check=FALSE)
stack3Excel(20,150, 2,2,2,2,2,2,9, existABC='no') #

#Note: To create an array, highlight the cells, then edit the formula, then type CTRL-SHIFT-ENTER.
#RApply("stack2",D6,D4, d9,d10,d11, f11, h6,h4)
#RApply("stack3",D22,D20, D25,D26,D27,d28,d29,d30,d31, f28,f29,f30,f31, h22,h20)
}
