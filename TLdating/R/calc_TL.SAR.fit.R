#' Estimation of the equivalent dose (ED) value for the SAR protocol
#' 
#' Internal function called by \link{analyse_TL.SAR}. \cr
#' This function estimates the equivalent dose (ED) based on a doses vector and a Lx/Tx vector provided. \cr
#' See details for more information.
#' 
#' @param doses
#'  \link{numeric} (\bold{required}): doses vector
#' @param LnTn
#'  \link{numeric} (\bold{required}): Ln/Tn.
#' @param LnTn.error
#'  \link{numeric} (\bold{required}): Error for the Ln/Tn.
#' @param LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx vector
#' @param LxTx.error
#'  \link{numeric} (\bold{required}): Error for the Lx/Tx vector
#' @param fitting.parameters
#'  \link{list} (with default): fitting parameters. See details.
#'
#' @details
#' This function estimates the equivalent dose based on the doses vector, Ln/Tn and the Lx/Tx matrix provided. \cr
#' Different fitting methods are available (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}). 
#' Moreover, the fitting can be weigthed or not. \cr
#' 
#' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.use.slope}}{
#'    \link{logical}: If the slope of the Q growth curve is reused for the sublinearity correction.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: Lowest regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: Highest regenerative dose used for the fitting.}
#' }
#' 
#' @return
#'  The function provides an \linkS4class{TLum.Results} object containing: \cr
#'  \describe{
#'    \item{\code{GC}}{
#'      \link{list}: fitting curve.}
#'    \item{\code{Q}}{
#'      \link{numeric}: equivalent dose estimation}
#'    \item{\code{Q.error}}{
#'      \link{numeric}: Error for the equivalent dose estimation}
#'    \item{\code{summary}}{
#'      \link{list}: parameters of the fitting result.}
#'  }
#' 
#' @seealso   
#'  \link{analyse_TL.SAR}.
#'  
#' @author David Strebler, University of Cologne (Germany).
#' 
## @export calc_TL.SAR.fit

calc_TL.SAR.fit <- function(
  
  doses,

  LnTn,
  
  LnTn.error,
    
  LxTx,
  
  LxTx.error,
  
  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE)
  
){  
  
  fit.method <- fitting.parameters$fit.method
  fit.weighted <- fitting.parameters$fit.weighted
  
  # Integrity Check ---------------------------------------------------------
  
  if (missing(LnTn)){
    stop("[calc_TL.SARfit] Error: LnTn is missing.")
  }
  
  if (missing(LnTn.error)){
    stop("[calc_TL.SARfit] Error: LnTn.error is missing.")
  }
  
  if (missing(LxTx)){
    stop("[calc_TL.SARfit] Error: LxTx is missing.")
  }
  
  if (missing(LxTx.error)){
    stop("[calc_TL.SARfit] Error: LxTx.error is missing.")
  }
  
  if (missing(doses)){
    stop("[calc_TL.SARfit] Error: doses is missing.")
  }
  # ------------------------------------------------------------------------------
    
  w <- 1/(LxTx.error^2)
  
  w[!is.finite(w)] <- 0.0001
  LxTx[!is.finite(LxTx)] <- 0.0001
  
  data <- data.frame(x=doses, 
                     y=LxTx)
  
  #LIN
  function.LIN <- function(a,b,x){a+b*x}
  
  #EXP
  function.EXP <- function(a,b,c,x){a*(1-exp(-(x+c)/b))}
  
  #EXP+LIN
  function.EXPLIN <- function(a,b,c,g,x){a*(1-exp(-(x+c)/b)+(g*x))}
  
  #EXP+EXP
  function.EXPEXP <- function(a1,a2,b1,b2,x){(a1*(1-exp(-(x)/b1)))+(a2*(1-exp(-(x)/b2)))}
  

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
    
    Q <- (LnTn-res["a"])/res["b"]
    
    temp.error <- sqrt(sum(LnTn.error^2,(res.error["a"])^2,na.rm = TRUE))
    temp.error.r <- temp.error/abs(LnTn-res["a"])
    Q.error.r <- sqrt(sum(temp.error.r^2,(res.error.r["b"])^2,na.rm=TRUE))
    Q.error <- abs(Q*Q.error.r)
    
  }else if(fit.method == "EXP"){
    
    # initial parameters
    # a
    a.i<-max(LxTx)
    
    # b
    temp.lm <- try(lm(log(data$y)~data$x),
                   silent = TRUE)
    
    if(class(temp.lm) =="try-error"){
      b.i<-1
    }else{
      b.i <- as.numeric(1/temp.lm$coefficients[2])
    } 
    
    #c
    temp.lm <- try(lm(data$y~data$x),
                   silent = TRUE)
    
    if(class(temp.lm) =="try-error"){
      c.i<-0
    }else{
      c.i <- as.numeric(abs(temp.lm$coefficients[1]/temp.lm$coefficients[2]))
    }  
      
    #Fitting
      # Weighted
    if(fit.weighted){       
      fit<-try(nls(y~function.EXP(a,b,c,x),
                   data=data,
                   start=c(a=a.i,
                           b=b.i,
                           c=c.i),
                   weights=w,
                   trace=FALSE,
                   algorithm="port",
                   nls.control(maxiter=500)),
               silent=TRUE)
      
      # Not weighted
    }else{
      fit<-try(nls(y~function.EXP(a,b,c,x),
                   data=data,
                   start=c(a=a.i,
                           b=b.i,
                           c=c.i),
                   trace=FALSE,
                   algorithm="port",
                   nls.control(maxiter=500)),
               silent=TRUE)
    }
    
    if (class(fit)=="try-error"){   
      
      print("pas ok")
      
      fit <- NA
      s <- list(a = NA,
                a.error = NA,
                b = NA,
                b.error = NA,
                c = NA,
                c.error = NA)
      Q <- 0
      Q.error <- NA
      
    }else{     
      print("ok")
      
      res <- summary(fit)$coefficients[,"Estimate"]
      res.error.r <- summary(fit)$coefficients[,"Std. Error"]
            
      res.error <- abs(res*res.error)
      
      s <- list(a = res["a"],
                a.error = res.error["a"],
                b = res["b"],
                b.error = res.error["b"],
                c = res["c"],
                c.error = res.error["c"])
      
      Q <- -res["c"]-res["b"]*log(1-LnTn/res["b"])
      Q.error <- 1 # temp
    }
    
  }else{
    stop("[calc_TL.SARfit] Error: fit.method not supported.")  
  } 

  result <- list(GC=fit,
                 Q=Q, 
                 Q.error=Q.error,
                 summary=s)
    
  new.TLum.Results.calc_TL.SAR.fit <- set_TLum.Results(data = result)

  return (new.TLum.Results.calc_TL.SAR.fit) 
}
  
