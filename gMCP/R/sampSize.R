#' Sample size calculations
#' 
#' Sample size calculations
#' 
#' @param graph A graph of class \code{\link{graphMCP}}.
#' @param esf ...
#' @param effSize ...
#' @param powerReqFunc One power requirement function or a list of these. 
#' % For example \code{function(x) {any(x)}} for at least one rejection or more complex user defined functions like \code{function(x) {(x[1]&&x[3])||(x[2]&&x[4])}}.
#' If one is interested in the power to reject hypotheses 1 and 3
#' one could specify: \cr\code{f=function(x) {x[1] && x[3]}}.\cr If the power
#' of rejecting hypotheses 1 and 2 is also of interest one would use a
#' (optionally named) list: \cr 
#' \code{f=list(power1and3=function(x) {x[1] && x[3]},}\cr
#' \code{power1and2=function(x) {x[1] && x[2]})}.
#' If the list has no names, the functions will be referenced 
#' to as "func1", "func2", etc. in the output.
#' @param target Target power that should be at least achieved. Either a numeric scalar between 0 and 1 or if parameter \code{powerReqFunc} is a list a numeric vector of the same length as \code{powerReqFunc}.
#' @param corr.sim Covariance matrix under the alternative.
#' @param corr.test Correlation matrix that should be used for the parametric test.
#' If \code{corr.test==NULL} the Bonferroni based test procedure is used. Can contain
#' NAs.
#' @param type What type of random numbers to use. \code{quasirandom} uses a
#' randomized Lattice rule, and should be more efficient than
#' \code{pseudorandom} that uses ordinary (pseudo) random numbers.
#' @param test In the parametric case there is more than one way to handle
#' subgraphs with less than the full alpha. If the parameter \code{test} is
#' missing, the tests are performed as described by Bretz et al. (2011), i.e.
#' tests of intersection null hypotheses always exhaust the full alpha level
#' even if the sum of weights is strictly smaller than one. If
#' \code{test="simple-parametric"} the tests are performed as defined in
#' Equation (3) of Bretz et al. (2011).
#' @param upscale Logical. If \code{upscale=FALSE} then for each intersection 
#' of hypotheses (i.e. each subgraph) a weighted test is performed at the 
#' possibly reduced level alpha of sum(w)*alpha, 
#' where sum(w) is the sum of all node weights in this subset.
#' If \code{upscale=TRUE} all weights are upscaled, so that sum(w)=1.
#' @param alpha ...
#' @param n.sim ...
#' @param verbose Logical, whether verbose output should be printed.
#' @param ... ...
#' @return ...
#' @examples
#' 
#' \dontrun{
#' graph <- BonferroniHolm(4)
#' powerReqFunc <- function(x) { (x[1] && x[2]) || x[3] }
#' #TODO Still causing errors / loops.
#' #sampSize(graph, alpha=0.05, powerReqFunc, target=0.8, mean=c(6,4,2) )
#' #sampSize(graph, alpha=0.05, powerReqFunc, target=0.8, mean=c(-1,-1,-1), nsim=100)
#' sampSize(graph, esf=c(1,1,1,1), effSize=c(1,1,1,1), 
#'          corr.sim=diag(4), powerReqFunc=powerReqFunc, target=0.8, alpha=0.05)
#' powerReqFunc=list('all(x[c(1,2)])'=function(x) {all(x[c(1,2)])},
#'                   'any(x[c(0,1)])'=function(x) {any(x[c(0,1)])})
#' sampSize(graph=graph, 
#'          effSize=list("Scenario 1"=c(2, 0.2, 0.2, 0.2), 
#'                       "Scenario 2"=c(0.2, 4, 0.2, 0.2)), 
#'          esf=c(0.5, 0.7071067811865476, 0.5, 0.7071067811865476),
#'          powerReqFunc=powerReqFunc, 
#'          corr.sim=diag(4), target=c(0.8, 0.8), alpha=0.025)
#'}
sampSize <- function(graph, esf, effSize, powerReqFunc, target,
                     corr.sim, alpha, corr.test = NULL,
                     type = c("quasirandom", "pseudorandom"),
                     upscale=FALSE, n.sim=10000, verbose=FALSE, ...) { # effSize, endPoints
  
  if (!is.list(effSize)) {
    effSize <- list(Scenario=effSize)
  }
  if (!is.list(powerReqFunc)) {
    powerReqFunc <- list(powerReqFunc=powerReqFunc)
  }
  # TODO Check whether names(...) has same length as ...
  
  # First determine 
  
  result <- list()
  df <- c()
  
  for (i in 1:length(powerReqFunc)) {
    for (scenario.name in names(effSize)) {
      prf.name <- names(powerReqFunc)[i]
      es <- effSize[[scenario.name]]
      prf <- powerReqFunc[[prf.name]]
      targFunc <- function(n, ...) {
        result <- calcPower(graph=graph, alpha=alpha, mean = es*sqrt(n)*esf,
                  corr.sim = corr.sim, corr.test = corr.test, type = type,
                  f=prf, upscale=upscale, ...)        
        return(result[[5]])
      }
      subResult <- sampSizeCore(100, targFunc=targFunc, alRatio=1, target=target[i], verbose=verbose, n.sim=n.sim)      
      df <- rbind(df, data.frame(Scenario=scenario.name, PowerFunc=prf.name, target=subResult$target, sampSize=subResult$samp.size))
      # n.sim = ceiling(n.sim/50)
      result[[length(result)+1]] <- subResult
      names(result)[[length(result)]] <- paste(prf.name, "-", scenario.name)
    }
  }
  return(df)
}



## Function for sample size calculation and functions to evaluate
## performance metrics for different sample sizes

#' Function for sample size calculation
#' 
#' Function for sample size calculation
#' 
#' For details see the manual and examples.
#' 
#' @param upperN \code{targFunc(upperN)} should be bigger than target (otherwise \code{upperN} is doubled until this is the case).
#' @param lowerN \code{targFunc(lowerN)} should be smaller than target (otherwise \code{lowerN} is halfed until this is the case).
#' @param targFunc The target (power) function that should be monotonically increasing in \code{n}.
#' @param target The target value. The function searches the \code{n} with \code{targFunc(n)-target<tol} and \code{targFunc(n)>target}.
#' @param tol Tolerance: The function searches the \code{n} with \code{targFunc(n)-target<tol} and \code{targFunc(n)>target}.
#' @param alRatio Allocation ratio.
#' @param Ntype Either \code{"arm"} or \code{"total"}.
#' @param verbose Logical, whether verbose output should be printed.
#' @param ... ...
#' @return Integer value \code{n} (of type numeric) with \code{targFunc(n)-target<tol} and \code{targFunc(n)>target}.
#' @author This function is taken from package DoseFinding under GPL from Bjoern Bornkamp, Jose Pinheiro and Frank Bretz
#' @examples 
#' 
#' f <- function(x){1/100*log(x)}
#' gMCP:::sampSizeCore(upperN=1000, targFunc=f, target=0.008, verbose=TRUE, alRatio=1)
#' 
sampSizeCore <- function (upperN, lowerN = floor(upperN/2),
                      targFunc, target, tol = 0.001, alRatio,
                      Ntype = c("arm", "total"), verbose = FALSE, ...){
  if (verbose) cat("Trying to find a sample size for power", target, "\n")
  
  ## target function to iterate
  func <- function(n, ...){
    targFunc(n, ...) - target
  }

  Ntype <- match.arg(Ntype)
  if (missing(alRatio)) stop("allocation ratios need to be specified")
  if (any(alRatio <= 0)) stop("all entries of alRatio need to be positive")
  
  alRatio <- alRatio/sum(alRatio)
  if(Ntype == "arm") {
    alRatio <- alRatio/min(alRatio)
  } 
  
  ## first call
  upper <- func(round(upperN*alRatio), ...)
  if(length(upper) > 1) stop(paste("targFunc(n) needs to evaluate to a vector of length 1, but returned:\n", paste(capture.output(dput(upper)), collapse="\n"), sep=""))
  if(!is.numeric(upper)) stop("targFunc(n) needs to evaluate to a numeric.")

  ## bracket solution
  # if (upper < 0) message("upper limit for sample size is raised")

  while (upper < 0) {
    upperN <- 2 * upperN
    upper <- func(round(upperN*alRatio), ...)
    message(paste("Upper limit for sample size is raised to ",upperN," (diff:",upper,")", sep=""))
  }
  
  lower <- func(round(lowerN*alRatio), ...)
  
  if (lower > 0) message("lower limit for sample size is decreased")

  while (lower > 0) {
    lowerN <- round(lowerN/2)
    if (lowerN == 0) stop("cannot find lower limit on n")
    lower <- func(round(lowerN*alRatio), ...)
  }

  ## now start bisection
  if (verbose) {
    cat("Upper N:", upperN, "Upper value", round(upper+target, 4), "\n")
    cat("Lower N:", lowerN, "Lower value", round(lower+target, 4), "\n\n")
  }
  
  current <- tol+1
  niter <- 0
  ## bisect sample size until tolerance is achieved
  while (abs(current) > tol & (upperN > lowerN + 1)) {
    currN <- round((upperN + lowerN)/2)
    current <- func(round(currN * alRatio), ...)
    if (current > 0) {
      upperN <- currN
    } else {
      lowerN <- currN
    }
    niter <- niter + 1
    if (verbose) {
      cat("Iter: ", niter, ", N = ", currN, ", current value = ",
          round(current+target, 4), "\n", sep = "")
    }
  }
  ## increase sample size so that the obtained value is larger than the target
  while (current < 0) {
    currN <- currN + 1
    current <- func(round(currN * alRatio), ...)
  }

  res <- list(samp.size = round(currN * alRatio),
              target = round(current+target, 4))
  attr(res, "alRatio") <- round(alRatio/min(alRatio), 4)
  attr(res, "target") <- target
  attr(res, "Ntype") <- Ntype
  class(res) <- "sampSize"
  res
}
