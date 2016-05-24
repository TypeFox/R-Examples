#' Estimation of the palaeodose (Q) with the MAAD protocol
#' 
#' Internal function called by \link{analyse_TL.MAAD}. \cr
#' This function makes a first estimation of the palaeodose based on a doses vector and a Lx/Tx vector provided. \cr
#' See details for more information.
#' 
#' 
#' 
#' @param LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx vector
#' @param LxTx.error
#'  \link{numeric} (\bold{required}): Error for the Lx/Tx vector
#' @param doses
#'  \link{numeric} (\bold{required}): doses vector
#' @param fitting.parameters
#'  \link{list} (with default): fitting parameters. See details.
#' 
#' @details
#' This function estimates the equivalent dose before any sublineary correction based on the doses vector and the Lx/Tx matrix provided. \cr
#' Different fitting methods are available (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}). 
#' Moreover, the fitting can be weigthed or not. \cr
#' 
#' #' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.use.slope}}{
#'    \link{logical}: If the slope of the Q growth curve is reused for the sublinearity correction.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: lowest regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: Highest regenerative dose used for the fitting.}
#' }
#' 
#' @return
#'  The function provides an \linkS4class{TLum.Results} object containing: \cr
#'  \describe{
#'    \item{\code{GC}}{
#'      \linkS4class{lm}: fitting result.}
#'    \item{\code{Q}}{
#'      \link{numeric}: equivalent dose estimation}
#'    \item{\code{Q.error}}{
#'      \link{numeric}: Error for the equivalent dose estimation}
#'    \item{\code{summary}}{
#'      \link{numeric}: parameters of the fitting result.}
#'  }
#'  
#'@seealso
#' \link{calc_TL.MAAD.fit.I},
#' \link{analyse_TL.MAAD}.
#' 
#' @author David Strebler, University of Cologne (Germany).
#' 
## @export calc_TL.MAAD.fit.Q

calc_TL.MAAD.fit.Q <- function(
  
  LxTx,
  
  LxTx.error,
  
  doses,
  
  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE)
  
){  
  
  V_METHODS <- c("LIN", "EXP","EXP+LIN", "EXP+EXP")
  
  #LIN
  function.LIN <- function(a,b,x){a+b*x}
  # a: horizontal translation
  # b: slope
  
  #EXP
  function.EXP <- function(a,b,h,x){a*(1-exp(-(x+h)/b))}
  # a: vectical scale, reflexion against x
  # b: horizontal scale, reflexion against y
  # h: horizontal translation
  
  #EXP+LIN
  function.EXPLIN <- function(a1,b1,h1,a2,x){a1*(1-exp(-(x+h1)/b1)+(a2*x))}
  #a1: (EXP comp) vectical scale, reflexion against x
  #b1: (EXP comp) horizontal scale, reflexion against y
  #h1: (EXP comp) horizontal translation
  #a2: (LIN comp) slope
  
  #EXP+EXP
  function.EXPEXP <- function(a1,b1,h1,a2,b2,h2,x){(a1*(1-exp(-(x+h1)/b1)))+(a2*(1-exp(-(x+h2)/b2)))}
  #a1: (EXP comp1) vectical scale, reflexion against x
  #b1: (EXP comp1) horizontal scale, reflexion against y
  #h1: (EXP comp1) horizontal translation
  #a2: (EXP comp2) vectical scale, reflexion against x
  #b2: (EXP comp2) horizontal scale, reflexion against y
  #h2: (EXP comp2) horizontal translation  
  
  
  # Integrity Check ---------------------------------------------------------
  
  if (missing(LxTx)){
    stop("[calc_TL.MAAD.fit] Error: LxTx is missing.")
  }
  
  if (missing(LxTx.error)){
    stop("[calc_TL.MAAD.fit] Error: LxTx.error is missing.")
  }
  
  if (missing(doses)){
    stop("[calc_TL.MAAD.fit] Error: doses is missing.")
  }
  # ------------------------------------------------------------------------------
  
  fit.method <- fitting.parameters$fit.method
  fit.weighted <- fitting.parameters$fit.weighted
  
  w <- 1/(LxTx.error^2)
  
  # ------------------------------------------------------------------------------
  # Check value
  if(is.null(fit.method)){
    stop("[calc_TL.MAAD.fit.I] Error: fit.method is missing")
  }else if(!(fit.method %in% V_METHODS)){
    stop("[calc_TL.MAAD.fit.I] Error: The Fitting method is not suported.")
  }
  
  if(!(is.logical(fit.weighted))){
    stop("[calc_TL.MAAD.fit.I] Error: fit.weighted is missing.")
  }
  
  if(length(LxTx) != length(LxTx.error)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx and LxTx.error have a different length.")
  }
  
  if(length(doses) != length(LxTx)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx and dose have a different length.")
  }
  
  if(FALSE %in% is.finite(LxTx)){
    LxTx[!is.finite(LxTx)] <- 0.0001
    #warning("[calc_TL.MAAD.fit.I] Warning: LxTx not always finite")
  }
  
  if(FALSE %in% is.finite(w)){
    w[!is.finite(w)] <- 0.0001   
    #warning("[calc_TL.MAAD.fit.I] Warning: w not always finite")
  }
  # ------------------------------------------------------------------------------
  
  data <- data.frame(x=doses, 
                     y=LxTx)
  
  if(fit.method == "LIN"){
    
    # Weighted
    if(fit.weighted){       
      fit <- lm(formula = y ~ x,
                data = data,
                weights=w)
      
      # Not weighted
    }else{
      fit <- lm(formula = y ~ x,
                data = data)
    }
    
    res <- summary(fit)$coefficients[,"Estimate"]
    names(res) <- c("a","b")

    res.error.r <- summary(fit)$coefficients[,"Std. Error"]
    names(res.error.r) <- c("a","b")

    res.error <- abs(res*res.error.r)
    
    s <- list(a = res["a"],
              a.error = res.error["a"],
              b = res["b"],
              b.error = res.error["b"])
    
    Q <- abs(res["a"]/res["b"])

    Q.error.r <- sqrt(sum(res.error.r^2,na.rm=TRUE))
    Q.error <- abs(Q.error.r*Q)
    
    
  }else if(fit.method == "EXP"){
    stop("[calc_TL.MAAD.fit] Error: EXP fitting not supported.")  
  
  }else{
    stop("[calc_TL.MAAD.fit] Error: fit.method not supported.")  
  } 
  
  result <- list(GC=fit,
                 Q=Q, 
                 Q.error=Q.error,
                 summary=s)
  
  new.TLum.Results.calc_TL.MAAD.fit <- set_TLum.Results(data = result)
  
  return (new.TLum.Results.calc_TL.MAAD.fit) 
}
