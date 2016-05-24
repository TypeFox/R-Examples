################################################# Class
setClass(  
  Class="ycInterExtra",
  representation=representation(
    matsin ="vector",
    observedvalues="vector",
    method="character",
    typeres ="character",
    coefficients="vector",
    matsout ="vector",
    fittedvalues="vector",
    residuals="vector",
    fwdrates ="vector",
    UFR="numeric",
    T_UFR="numeric",
    extrapvalues="vector"
  )
)

################################################## Getter - Setter

# maturities
setGeneric(
  name="setMatsIn",
  def = function(.Object, x)
  {
    standardGeneric("setMatsIn")
  }
)
setMethod (f = "setMatsIn",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, x)
             {
               .Object@matsin <- x
               return(.Object)
             })
setGeneric(
  name="getMatsIn",
  def = function(.Object)
  {
    standardGeneric("getMatsIn")
  }
)
setMethod (f = "getMatsIn",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@matsin)
             })


setGeneric(
  name="setMatsOut",
  def = function(.Object, x)
  {
    standardGeneric("setMatsOut")
  }
)
setMethod (f = "setMatsOut",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, x)
             {
               .Object@matsout <- x
               return(.Object)
             })

setGeneric(
  name="getMatsOut",
  def = function(.Object)
  {
    standardGeneric("getMatsOut")
  }
)
setMethod (f = "getMatsOut",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@matsout)
             })


# observed values
setGeneric(
  name="setObservedvalues",
  def = function(.Object, x)
  {
    standardGeneric("setObservedvalues")
  }
)
setMethod (f = "setObservedvalues",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, x)
             {
               .Object@observedvalues <- x
               return(.Object)
             })

setGeneric(
  name="getObservedvalues",
  def = function(.Object)
  {
    standardGeneric("getObservedvalues")
  }
)
setMethod (f = "getObservedvalues",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@observedvalues)
             })

# coefficients
setGeneric(
  name="setCoefficients",
  def = function(.Object, x)
  {
    standardGeneric("setCoefficients")
  }
)
setMethod (f = "setCoefficients",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, x)
             {
               .Object@coefficients <- unlist(x)
               return(.Object)
             })
setGeneric(
  name="getCoefficients",
  def = function(.Object)
  {
    standardGeneric("getCoefficients")
  }
)
setMethod (f = "getCoefficients",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@coefficients)
             })

# residuals
setGeneric(
  name="setResiduals",
  def = function(.Object, x)
  {
    standardGeneric("setResiduals")
  }
)
setMethod(f = "setResiduals",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, x)
             {
               .Object@residuals <- x
               return(.Object)
             })
setGeneric(
  name="getResiduals",
  def = function(.Object)
  {
    standardGeneric("getResiduals")
  }
)
setMethod (f = "getResiduals",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@residuals)
             })

# Forward rates
setGeneric(
  name="setFwdrates",
  def = function(.Object, x)
  {
    standardGeneric("setFwdrates")
  }
)
setMethod (f = "setFwdrates",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, x)
             {
               .Object@fwdrates <- x
               return(.Object)
             })
setGeneric(
  name="getFwdrates",
  def = function(.Object)
  {
    standardGeneric("getFwdrates")
  }
)
setMethod (f = "getFwdrates",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@fwdrates)
             })


# fitted values
setGeneric(
  name="setFittedvalues",
  def = function(.Object, yM = NULL, p = NULL, matsin, matsout, 
                 method=c("NS", "SV", "SW", "HCSPL"), typeres=c("rates", "prices"))
  {
    standardGeneric("setFittedvalues")
  }
)  
              
setMethod (f = "setFittedvalues",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, yM = NULL, p = NULL, matsin, matsout, 
                      method=c("NS", "SV", "SW", "HCSPL"), typeres=c("rates", "prices"))
             {
               .Object@matsin <- matsin               
               .Object@matsout <- matsout
               method <- match.arg(method)
               typeres <- match.arg(typeres)
               .Object@typeres <- typeres               
               .Object@method <- method
               J <- length(matsout)
               indicemat <- pmatch(matsin, matsout)
               
               if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))                 
               {
                     .Object <- setObservedvalues(.Object, p)
                     res <- ycInterpolation(p = p, matsin = matsin, matsout = matsout, 
                                            method=method, typeres=typeres)
                     .Object <- setCoefficients(.Object, res$coefficients)                 
                     .Object@fittedvalues <- res$values
                     .Object <- setFwdrates(.Object, res$fwd)
                     
                     if (typeres == "prices")
                     { 
                       x <- p - res$values[indicemat]
                     }
                     else 
                    {
                            P <- pricefromeuribor(0, matsout, res$values)
                            x <- p - P[indicemat]
                     }
                     
                     .Object <- setResiduals(.Object, x)           
               }
               
               if ((is.null(p) || missing(p)) && (!is.null(yM) || !missing(yM)))                 
               {
                   .Object <- setObservedvalues(.Object, yM)                   
                   res <- ycInterpolation(yM = yM, matsin = matsin, matsout = matsout, 
                                          method=method, typeres=typeres)                 
                   .Object <- setCoefficients(.Object, res$coefficients)                 
                   .Object@fittedvalues <- res$values
                   .Object <- setFwdrates(.Object, res$fwd)
                   
                     if (typeres == "rates")
                     {
                       P <- pricefromeuribor(0, matsout, res$values)
                       x <- yM - res$values[indicemat]
                     }
                     else {                       
                       x <- yM - euriborfromprice(0, matsout[indicemat], res$values[indicemat])
                     }
                     
                   .Object <- setResiduals(.Object, x)
               }
               
               return(.Object)
             })

setGeneric(
  name="getFittedvalues",
  def = function(.Object)
  {
    standardGeneric("getFittedvalues")
  }
)
setMethod (f = "getFittedvalues",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@fittedvalues)
             })

# extrapolated values
setGeneric(
  name="setExtrapvalues",
  def = function(.Object, yM = NULL, p = NULL, matsin, matsout, 
                 method=c("NS", "SV", "SW"), typeres=c("rates", "prices"), UFR, 
                 T_UFR = NULL)
  {
    standardGeneric("setExtrapvalues")
  }
)  

setMethod (f = "setExtrapvalues",
           signature = "ycInterExtra", 
           definition = 
             function(.Object, yM = NULL, p = NULL, matsin, matsout, 
                      method=c("NS", "SV", "SW"), typeres=c("rates", "prices"), UFR, 
                      T_UFR)
             {
               .Object <- setMatsIn(.Object, matsin)               
               .Object <- setMatsOut(.Object, matsout)                
               method <- match.arg(method)
               typeres <- match.arg(typeres)
               .Object@method <- method
               .Object@typeres <- typeres      
               .Object@UFR <- UFR               
               indicemat <- pmatch(matsin, matsout)
               J <- length(matsout)
                              
               if(!is.null(T_UFR)) 
               {
                 if ((method != "SW")) warning("unused parameter T_UFR")
                 .Object@T_UFR <- T_UFR
               }
                              
               if (!is.null(p) && is.null(yM))                 
               {
                 .Object <- setObservedvalues(.Object, p)   
                  res <- ycExtrapolation(p = p, matsin = matsin, matsout = matsout, 
                                          method = method, typeres = typeres, UFR = UFR, T_UFR = T_UFR)                 
                 .Object@fittedvalues <- res$values
                 .Object <- setFwdrates(.Object, res$fwd)                                  
                 
                 if (typeres == "prices")
                 {
                   x <- p - res$values[indicemat]
              }
                 else 
                  {
                   P <- pricefromeuribor(0, matsout, res$values)
                   x <- p - P[indicemat]
                }
                                    
                 .Object <- setResiduals(.Object, x)                 
               }
               
               if (is.null(p) && !is.null(yM))                 
               {
                 .Object <- setObservedvalues(.Object, yM)
                 
                 res <- ycExtrapolation(yM = yM, matsin = matsin, matsout = matsout, 
                                          method = method, typeres = typeres, UFR = UFR, T_UFR = T_UFR)
                 .Object@fittedvalues <- res$values
                 .Object <- setFwdrates(.Object, res$fwd)
                 
                 if (typeres == "rates")
                 {
                   x <- yM - res$values[indicemat]
                   P <- pricefromeuribor(0, matsout, res$values)
              }
                 else {
                   x <- yM - euriborfromprice(0, matsout[indicemat], res$values[indicemat])
                 }                 
                 .Object <- setResiduals(.Object, x)                 
               }
                              
               .Object <- setCoefficients(.Object, res$coefficients)
               
               .Object@extrapvalues <- res$values               
               
               return(.Object)
             })

setGeneric(
  name="getExtrapvalues",
  def = function(.Object)
  {
    standardGeneric("getExtrapvalues")
  }
)
setMethod (f = "getExtrapvalues",
           signature = "ycInterExtra", 
           definition = 
             function(.Object)
             {
               return(.Object@extrapvalues)
             })

########################################################## Fonctions de l'interface

ycinter <- function(yM = NULL, p = NULL, matsin, matsout, 
                    method=c("NS", "SV", "SW", "HCSPL"), 
                    typeres=c("rates", "prices"))
{
    if ((missing(yM) || is.null(yM)) && (missing(p) || is.null(p))) stop("no input values")
    
    y <- new("ycInterExtra")      
    method <- match.arg(method)
    typeres <- match.arg(typeres)
    
    if (!missing(yM) || !is.null(yM))
    {
      y <- setFittedvalues(y, yM = yM, matsin = matsin, matsout = matsout, 
                         method=method, typeres=typeres)
    }
    
    if (!missing(p) || !is.null(p))
    {
      y <- setFittedvalues(y, p = p, matsin = matsin, matsout = matsout, 
                           method=method, typeres=typeres)
    }
        
    return(y)
}

ycextra <- function(yM = NULL, p = NULL, matsin, matsout, 
                    method=c("NS", "SV", "SW"), 
                    typeres=c("rates", "prices"), UFR, 
                    T_UFR = NULL)
{
  if ((missing(yM) || is.null(yM)) && (missing(p) || is.null(p))) stop("no input values")
  
  y <- new("ycInterExtra")      
  method <- match.arg(method)
  typeres <- match.arg(typeres)
  
  if (!missing(yM) || !is.null(yM))
  {
    y <- setExtrapvalues(y, yM = yM, matsin = matsin, matsout = matsout, 
                         method=method, typeres=typeres, UFR = UFR, T_UFR = T_UFR)
  }
  
  if (!missing(p) || !is.null(p))
  {
    y <- setExtrapvalues(y, p = p, matsin = matsin, matsout = matsout, 
                          method=method, typeres=typeres, UFR = UFR, T_UFR = T_UFR)
  }
  
  return(y)
}

 coeffs <- function(.Object)
 {
   return(getCoefficients(.Object))
 }

 deviance <- function(.Object)
 {
   return(crossprod(getResiduals(.Object))[1,1])
 }

 residuals<- function(.Object)
 {
   u <- .Object@matsin
   return(ts(getResiduals(.Object), deltat=u[2]-u[1]))
 }

fitted<- function(.Object) 
{
  t <- getMatsOut(.Object)
  return(ts(getFittedvalues(.Object), deltat=t[2]-t[1]))
}

ycsummary <-   function(.Object)
{
   y <- list(obs = .Object@observedvalues, matsin = .Object@matsin, 
   coeff = .Object@coefficients, fitted =.Object@fittedvalues, 
   matsout = .Object@matsout, resid = .Object@residuals, 
   typeres = .Object@typeres)
   
   if(max(y$obs) < 0.2 && (y$typeres != "rates"))
   {
     y$fitted <- euriborfromprice(0, y$matsout, y$fitted) 
   }
     
   if(max(y$obs) > 0.2 && (y$typeres != "prices"))
   {
     y$fitted <- pricefromeuribor(0, y$matsout, y$fitted)
   }       
               
   n <- length(y$obs)
   if (.Object@method == "SW")
   {p <- length(y$coeff) - 1}
   else {p <- length(y$coeff)}
   
   cat("Residuals:", "\n")
   print(summary(y$resid))
   cat("\n")
               
   cat("Coefficients:", "\n")
   cat(y$coeff, "\n")
   cat("\n")
   
   # Total sum of squares
   yobs <- y$obs
   ybar <- mean(y$obs)
   typeres <- y$typeres
   u <- y$matsin
   t <- y$matsout
   indicemat <- pmatch(u, t)
   
   ychapeau <- y$fitted[indicemat]
   
   resid <- y$resid
   
   ESS <- crossprod(rep.int(ybar, n) - ychapeau)[1,1]
   SSR <- crossprod(yobs - ychapeau)[1,1]
   SST <- ESS + SSR
   
    cat("Total sum of squares:", "\n")
    cat(SST, "\n")                             
    cat("with", n - 1,"degrees of freedom","\n")      
    cat("\n")
               
    cat("Explained sum of squares:", "\n")
    cat(ESS, "\n")               
    cat("with", p - 1,"degrees of freedom","\n")      
    cat("\n")
               
    cat("Residual sum of squares:", "\n")
    cat(SSR, "\n")
    cat("with", n - p,"degrees of freedom","\n")
    cat("\n")
               
    Rsquared <- 1 - SSR/SST
    cat("Multiple R-squared ", "*", "Adjusted R-squared", "\n")
    cat(Rsquared, "*", 1 - ((n-1)/(n-p))*(1 - Rsquared), "\n")
    cat("\n")               
               
    if (try(shapiro.test(y$resid)$p.value, TRUE) >= 0.05 || n >= 30)
    {
     Fstat <- ((n-p)/(p-1))*(Rsquared/(1 - Rsquared))
     suppressWarnings(cat("F-statistic:", Fstat, "on", 
     n - p,"and", p-1, "degrees of freedom, p-value:", pf(Fstat, n-p, p-1), "\n"))
    }
}

ycplot <- function(.Object){               
               
indicatrice <- pmatch(.Object@matsin, .Object@matsout)               

par(mfrow=c(2,2))

ord <- .Object@fittedvalues
x <- .Object@observedvalues

    if(.Object@typeres == "rates")
    {       
                      if (max(.Object@observedvalues) > 0.2) 
                      {x <- euriborfromprice(0, .Object@matsin, .Object@observedvalues)}                  
                       
                       plot(x = .Object@matsout, y = ord, 
                            main = "(Red) Observed values and (Black) Fitted values", xlab = "Maturity",
                            ylab = "Yield to maturity")               
                       points(x = .Object@matsin, y = x, col="red", pch=3)
                        
                       plot(x = x, y = ord[indicatrice], 
                            main = "Observed vs Fitted values", xlab="observed values", ylab="fitted values")
                       abline(0, 1, col="red")                                    
    } 
    else 
    { 
                     if (max(.Object@observedvalues) < 0.2) 
                     {x <- pricefromeuribor(0, .Object@matsin, .Object@observedvalues)}                  
                     
                       plot(x = .Object@matsout, y = ord, 
                            main = "(Red) Observed values and (Black) Fitted values", xlab = "Maturity",
                            ylab = "Zero-coupon price")               
                       points(x = .Object@matsin, y = x, col="red", pch=3)
                     
                       plot(x = x, y = ord[indicatrice], 
                            main = "Observed vs Fitted values", xlab="observed values", ylab="fitted values")
                       abline(0, 1, col="red")                 
    }               
               
               hist(.Object@residuals, main = paste("Histogram and density", "\n","of residuals"), 
                    xlab = "residuals")
               lines(density(.Object@residuals), col="red")

              if(max(.Object@matsout) <= max(.Object@matsin))
              {
               qqnorm(.Object@residuals, main = "Residuals' Normal Q-Q plot")
               qqline(.Object@residuals, col="red")
              }
              else
              {
                plot(forwardrates(.Object), main = "Extrapolated forward rates")
              }
}


forwardrates <- function(.Object)
{
  t <- .Object@matsout
  return(ts(getFwdrates(.Object), deltat=t[2]-t[1]))
}

as.list <- function(.Object)
{
   y <- list(matsin = .Object@matsin, obs = .Object@observedvalues, method = .Object@method, typeres = .Object@typeres,
             coeff = .Object@coefficients, matsout = .Object@matsout, fitted =.Object@fittedvalues, fwdrates = .Object@fwdrates,
             resid = .Object@residuals, UFR = .Object@UFR, T_UFR = .Object@T_UFR, extra = .Object@extrapvalues)
   return(y)
}