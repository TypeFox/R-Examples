#' Estimation of the sublinearity value for the MAAD protocol
#' 
#' Internal function called by \link{analyse_TL.MAAD}. \cr
#' This function estimates the sublinearity correction based on the dose vector and the Lx/Tx vector provided. \cr
#' See details for more information.
#' 
#' @param LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx vector
#' @param LxTx.error
#'  \link{numeric} (\bold{required}): Error for the Lx/Tx vector
#' @param doses
#'  \link{numeric} (\bold{required}): doses vector
#' @param slope
#'  \link{list} (with default): Property of the additive growth curve.
#' @param fitting.parameters
#'  \link{list} (with default): fitting parameters. See details.
#' 
#' @details
#' 
#' This function estimates the sublinearity correction based on the doses vector and the Lx/Tx matrix provided. \cr
#' Different fitting methods are available (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}). 
#' Morover, the fitting can be weigthed or not. \cr
#' If the fitting parameter \code{fit.use.slope} is \code{TRUE}, the function will use the data 
#' from \code{slope} to define the fitting curve for the sublinearity correction.
#' In that case, the sublinearity correction growth curve will be parallel to the additive growth curve.    
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
#'      \linkS4class{lm}: The fitting result.}
#'    \item{\code{i}}{
#'      \link{numeric}: The sublinearity correction estimation for the given equivalent dose}
#'    \item{\code{I.error}}{
#'      \link{numeric}: The error for the sublinearity correction estimation}
#'    \item{\code{summary}}{
#'      \link{numeric}: The parameters of the fitting result.}
#'  }
#'  
#' @seealso
#'  \link{calc_TL.MAAD.fit.Q},
#'  \link{analyse_TL.MAAD}.
#'  
#' @author David Strebler, University of Cologne (Germany).
#' 
## @export calc_TL.MAAD.fit.I

calc_TL.MAAD.fit.I <- function(
  
  LxTx,
  
  LxTx.error,
  
  doses,
  
  slope = NULL,
  
  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE,
                          fit.use.slope=FALSE)
    
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
  function.EXPEXP <- function(a1,a2,h1,b1,b2,h2,x){(a1*(1-exp(-(x+h1)/b1)))+(a2*(1-exp(-(x+h2)/b2)))}
  #a1: (EXP comp1) vectical scale, reflexion against x
  #b1: (EXP comp1) horizontal scale, reflexion against y
  #h1: (EXP comp1) horizontal translation
  #a2: (EXP comp2) vectical scale, reflexion against x
  #b2: (EXP comp2) horizontal scale, reflexion against y
  #h2: (EXP comp2) horizontal translation  
  
  # Integrity Check ---------------------------------------------------------  
  if (missing(doses)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx is missing.")
  }else if(!is.numeric(LxTx)){
    stop("[calc_TL.MAAD.fit.I] Error: doses is not of type 'numeric'.")
  }
  
  if (missing(LxTx.error)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx.error is missing.")
  }else if(!is.numeric(LxTx.error)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx.error is not of type 'numeric'.")
  }
  
  if (missing(doses)){
    stop("[calc_TL.MAAD.fit.I] Error: doses is missing.")
  }else if(!is.numeric(doses)){
    stop("[calc_TL.MAAD.fit.I] Error: doses is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------
  
  fit.method <- fitting.parameters$fit.method
  fit.weighted <- fitting.parameters$fit.weighted
  fit.use.slope <- fitting.parameters$fit.use.slope

  w <- 1/(LxTx.error^2)
  
  # ------------------------------------------------------------------------
  # Check value
  if(is.null(fit.method)){
    stop("[calc_TL.MAAD.fit.I] Error: fit.method is missing")
  }else if(!(fit.method %in% V_METHODS)){
    stop("[calc_TL.MAAD.fit.I] Error: The Fitting method is not suported.")
  }
  
  if(!(is.logical(fit.weighted))){
    stop("[calc_TL.MAAD.fit.I] Error: fit.weighted is missing.")
  }

  if(!(is.logical(fit.use.slope))){
    stop("[calc_TL.MAAD.fit.I] Error: fit.use.slope is missing.")
  }
  
  if(length(LxTx) != length(LxTx.error)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx and LxTx.error have a different length.")
  }

  if(length(doses) != length(LxTx)){
    stop("[calc_TL.MAAD.fit.I] Error: LxTx and dose have a different length.")
  }
  
  if(fit.use.slope){
    if(is.null(slope)){
      fit.use.slope <- FALSE
      warning("[calc_TL.MAAD.fit.I] Warning: slope is missing.")
    }
  }
  
  if(FALSE %in% is.finite(LxTx)){
    LxTx[!is.finite(LxTx)] <- 0.0001
    #warning("[calc_TL.MAAD.fit.I] Warning: LxTx not always finite")
  }
  
  if(FALSE %in% is.finite(w)){
    w[!is.finite(w)] <- 0.0001   
    #warning("[calc_TL.MAAD.fit.I] Warning: w not always finite")
  }
  # ------------------------------------------------------------------------

  error.r <- 0.1
  if(fit.use.slope){
        
    if(fit.method == "LIN"){
      b <- slope$b
      b.error <- slope$b.error
      
    }else if(fit.method == "EXP"){
      a <- slope$a
      b <- slope$b
      
      a.error <- abs(a*error.r)
      b.error <- abs(b*error.r)
      
    }else if(fit.method == "EXP+LIN"){
      a1 <- slope$a1
      b1 <- slope$b1
      a2 <- slope$a2
      
      a1.error <- abs(a1*error.r)
      b1.error <- abs(b1*error.r)
      a2.error <- abs(a2*error.r)
      
    }else if(fit.method == "EXP+EXP"){
      a1 <- slope$a1
      b1 <- slope$b1
      a2 <- slope$a2
      b2 <- slope$b2
      
      a1.error <- abs(a1*error.r)
      b1.error <- abs(b1*error.r)
      a2.error <- abs(a2*error.r)
      b2.error <- abs(b2*error.r)  
    }
  }

  
  # ------------------------------------------------------------------------
  # Check value
  if(fit.use.slope){
    if(fit.method == "LIN"){
      if(!is.finite(b)){
        stop("[calc_TL.MAAD.fit.I] Error: b is not finite.")
      }
      
    }else if(fit.method == "EXP"){
      if(!is.finite(a)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: a is not finite.")
      }else if(!is.finite(b)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: b is not finite.")
      }
      
    }else if(fit.method == "LIN+EXP"){
      if(!is.finite(a1)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: a1 is not finite.")
      }else if(!is.finite(b1)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: b1 is not finite.")
      }else if(!is.finite(a2)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: a2 is not finite.")
      }
      
    }else if(fit.method == "EXP+EXP"){
      if(!is.finite(a1)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: a1 is not finite.")
      }else if(!is.finite(b1)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: b1 is not finite.")
      }else if(!is.finite(a2)){
        warning("[calc_TL.MAAD.fit.I] Error: a2 is not finite.")
      }else if(!is.finite(b2)){
        fit.use.slope <- FALSE
        warning("[calc_TL.MAAD.fit.I] Error: b2 is not finite.")
      }
    }
  }
  # ------------------------------------------------------------------------
  
  data <- data.frame(x=doses, 
                     y=LxTx)
  
  if(fit.method == "LIN"){
    
    if(fit.use.slope){      
      temp.a <- LxTx - b*doses
      temp.a.error <- sqrt(LxTx.error^2 + (b.error*doses)^2)
      
      temp.a.w <- 1/(temp.a.error^2)
      
      if(fit.weighted){
        # Weighted
        a <- sum(temp.a.w*temp.a,na.rm = TRUE)/sum(temp.a.w)
        a.error <- 1/sqrt(sum(temp.a.w))
        
      }else{
        # Not weighted
        a <- mean(temp.a)
        a.error <- mean(temp.a.error)
      }
      
      #fit
      fit <- c(a,b)
      names(fit) <- c("a","b")
      
    }else{
      # No slope
      if(fit.weighted){  
        # Weighted
        fit <- lm(formula = y ~ x,
                  data = data,
                  weights=w)
        
      }else{
        # Not weighted
        fit <- lm(formula = y ~ x,
                  data = data)
      }
      
      res <- summary(fit)$coefficients[,"Estimate"]
      res.error.r <- summary(fit)$coefficients[,"Std. Error"]
      
      a <- res[1]
      a.error <- res.error.r[1]*a
      b <- res[2]
      b.error <- res.error.r[2]*b
    }
    
    s <- list(a = a,
              a.error = a.error,
              b = b,
              b.error = b.error)
    
    I <- abs(a/b)
    
    I.error.r <- sqrt((a.error/a)^2+(b.error/b)^2)
    I.error <- abs(I.error.r*I)
    
    
  }else if(fit.method == "EXP"){
    stop("[calc_TL.MAAD.fit.I] Error: EXP fitting not supported yet.")  
    
  }else if(fit.method == "EXP+LIN"){
    stop("[calc_TL.MAAD.fit.I] Error: EXP+LIN fitting not supported yet.")  
  
  }else if(fit.method == "EXP+EXP"){
    stop("[calc_TL.MAAD.fit.I] Error: EXP+EXP fitting not supported yet.")  
  
  }else{
    stop("[calc_TL.MAAD.fit.I] Error: fit.method not supported.")  
  } 
  
  if(!is.finite(I)){
    I <- NA
    I.error <- NA
    warning("[calc_TL.MAAD.fit.I] Warning: I is not finite.")  
  }
  
  result <- list(GC=fit,
                 I=I, 
                 I.error=I.error,
                 summary=s)
  
  new.TLum.Results.calc_TL.MAAD.fit.I <- set_TLum.Results(data = result)
  
  return (new.TLum.Results.calc_TL.MAAD.fit.I) 
}
