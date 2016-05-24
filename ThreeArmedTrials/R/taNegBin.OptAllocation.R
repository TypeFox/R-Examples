#' @title Optimal sample size for three-armed clinical trials with negative binomial distributed endpoints
#' @description Calculate optimal sample size allocation for Wald-type test for superiority or non-inferiority of the experimental treatment versus reference treatment with respect to placebo
#' @param rateExp A numeric value specifying the rate of the experimental treatment group in the alternative hypothesis
#' @param rateRef A numeric value specifying the rate of the reference treatment group in the alternative hypothesis
#' @param ratePla A numeric value specifying the rate of the placebo treatment group in the alternative hypothesis
#' @param shape A numeric value specifying the shape parameter
#' @param Delta A numeric value specifying the non-inferiority/superiority margin
#' @param type A character string determing how the variance for the Wald-type test statistic is estimated, must be \emph{restricted}, or \emph{unrestricted} 
#' @param n The total sample size. This parameter is only mandatory for \emph{type='restricted'}. For \emph{type='unrestricted'}, this parameter is optional.
#' @param sig.level A numeric value specifying the significance level (type I error probability). This parameter is only mandatory for \emph{type='restricted'}. For \emph{type='unrestricted'}, this parameter is optional.
#' @return A list with class "power.htest" containing the following components:
#' \item{n}{The total sample size. Not mandatory.}
#' \item{Delta}{A numeric value specifying the non-inferiority/superiority margin}
#' \item{type}{A character string indicating what type of Wald-type test will be performed}
#' \item{allocation}{A vector with the sample size allocation (nExp/n, nRef/n, nPla/n)}
#' \item{rateExp}{A numeric value specifying the rate of the experimental treatment group in the alternative hypothesis}
#' \item{rateRef}{A numeric value specifying the rate of the reference treatment group in the alternative hypothesis}
#' \item{ratePla}{A numeric value specifying the rate of the placebo treatment group in the alternative hypothesis}
#' \item{shape}{A numeric value specifying the shape parameter}
#' \item{nExp}{A numeric value specifying the number of sample in the experimental treatment group}
#' \item{nRef}{A numeric value specifying the number of sample in the reference treatment group}
#' \item{nPla}{A numeric value specifying the number of sample in the placebo treatment group}
#' @examples 
#' # Example for type = 'unrestricted' 
#' taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  type = 'unrestricted', n = 1048, sig.level = 0.025)
#' taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  type = 'unrestricted')
#' 
#' # Example for type = 'restricted'. 
#' \dontrun{
#' taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, 
#'  type = 'restricted', n = 500, sig.level = 0.025)
#' }
#' @export
#' @keywords optimalSampleSizeAllocation NegativeBinomial
taNegBin.OptAllocation <- function(rateExp, rateRef, ratePla, shape, Delta, type = c('restricted', 'unrestricted'), n = NULL, sig.level = NULL){

  type <- match.arg(type)
  if( missing(rateExp) ){stop("'rateExp' is missing.")}
  if( missing(rateRef) ){stop("'rateRef' is missing.")}
  if( missing(ratePla) ){stop("'ratePla' is missing.")}
  if( missing(shape) ){stop("'shape' is missing.")}
  if( missing(Delta) ){stop("'Delta' is missing.")}
  
  # Calculate effect size and check if parameters are located in the alternative
  effect <- rateExp - Delta * rateRef + (Delta - 1) * ratePla  
  if( effect >= 0 ){
    stop('Parameter vector is not located in the alternative.')
  }

  
  sigmaExp <- sqrt( rateExp * (1 + rateExp * shape) )
  sigmaRef <- sqrt( rateRef * (1 + rateRef * shape) )
  sigmaPla <- sqrt( ratePla * (1 + ratePla * shape) )
  
  if( type == 'unrestricted' ){
    w <- c(sigmaExp, Delta * sigmaRef, abs(1-Delta)*sigmaPla) / (sigmaExp + Delta * sigmaRef + abs(1-Delta)*sigmaPla)
    method <- 'Optimal sample size allocation for Wald-type test with unrestricted variance estimation under assumption of negatively binomial distributed endpoints.'

    if( !is.null(n) && is.numeric(n) && (n>=7)){
      nExp <- round(w[1] * n)
      nRef <- round(w[2] * n)
      nPla <- n - nExp - nRef
      w <- c(nExp, nRef, nPla)/n
      if( !is.null(sig.level) && (is.numeric(sig.level) && (sig.level>0) && (sig.level<1)) ){
        VarUnres <- sigmaExp^2 / w[1] + Delta^2 * sigmaRef^2 / w[2] + (1-Delta)^2 * sigmaPla^2 / w[3]
        power <- pnorm(qnorm(sig.level) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1)
      }
      
      return(
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
                       nPla = nPla),
                  class = "power.htest")
      )
    }
    else{
      return(
        structure(list(method = method,
                       rateExp = rateExp,
                       rateRef = rateRef,
                       ratePla = ratePla,
                       shape = shape,
                       type = type,
                       Delta = Delta,
                       allocation = w),
                  class = "power.htest")
      )
    }
    
  }


     

  if( type == 'restricted' ){

    if( is.null(n) || !is.numeric(n) ){
      stop("'n' is missing.")
    }
    if( n < 7 ){
      stop("'n' has to be larger than 6.")
    }
    
    if( is.null(sig.level) || !is.numeric(sig.level) ){
      stop("'sig.level' is missing.")
    }
    if( (sig.level<=0) || (sig.level>=1) ){
      stop("'sig.level' has to be between 0 and 1.")
    }
    
    # Function to be minimized for optimal sample size
    PowerRmlAlloc <- function(w){
      wIn <- w
      w <- c(wIn, 1-wIn[1]-wIn[2])
      VarRes <- taNegbin.LimitRestMLE(rateExp1=rateExp, rateRef1=rateRef, 
                            ratePla1=ratePla, shape1=shape, Delta=Delta, 
                            allocation = w
      )$sigma2.rest
      VarUnres <- sigmaExp^2 / w[1] + Delta^2 * sigmaRef^2 / w[2] + (1-Delta)^2 * sigmaPla^2 / w[3]
      Power <- pnorm(qnorm(sig.level)*sqrt(VarRes)/ sqrt(VarUnres) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1)
      return(-Power)
    }
    
    # Define restrictions and calculate optimal allocation for unrestricted case as starting value
    ui <- matrix(c(c(1,1),diag(1, nrow=2, ncol = 2)), byrow = T, ncol = 2)
    wUnres <- c(sigmaExp, Delta * sigmaRef, abs(1-Delta)*sigmaPla) / (sigmaExp + Delta * sigmaRef + abs(1-Delta)*sigmaPla)
    # Caculation optimal allocation
    wOptRML <- constrOptim( theta=wUnres[1:2] , f=PowerRmlAlloc, grad = NULL,ui = ui, ci=c(-1,0,0))$par

    # Grid search around theoretical optimal allocation for feasiable optimal allocation

    nExp <- round(wOptRML[1]*n)
    nRef <- round(wOptRML[2]*n)
    nPla <- n - nExp - nRef
    w <- c(nExp, nRef, nPla) / n
    VarUnres <- sigmaExp^2 / w[1] + Delta^2 * sigmaRef^2 / w[2] + (1-Delta)^2 * sigmaPla^2 / w[3]
    VarRes <- taNegbin.LimitRestMLE(rateExp1=rateExp, rateRef1=rateRef, 
                                    ratePla1=ratePla, shape1=shape, Delta=Delta, 
                                    allocation = w
    )$sigma2.rest
    power <- pnorm(qnorm(sig.level)*sqrt(VarRes)/ sqrt(VarUnres) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1)
    
    # Define 'method' for output
    method <- 'Optimal sample size allocation for Wald-type test with maximum-likelihood variance estimation under assumption of negatively binomial distributed endpoints.'
    
    return(
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
                     nPla = nPla),
                class = "power.htest")
    )
    }
}
