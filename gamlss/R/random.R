#-------------------------------------------------------------------------------
# this function is a modifyied version of  Hastie's S-plus gam test function
# it allows a random effect fit for a factor
# TO DO:  i) check for weighted out obseervations
#        ii) put prediction
#       iii) check if weighted levels 
random <- function(x, df = NULL, lambda = NULL, start=10) 
{
  scall <- deparse(sys.call(), width.cutoff = 500L)
 if(!inherits(x, "factor")) # | !is.category(xvar))
    stop("random() expects a factor as its first argument")
 if (!is.null(df))
   {
     nlevelsM <- nlevels(x)
     df <- if(df>=nlevelsM) nlevelsM  else df
     df <- if(df<=0) 0.001  else df
   }
   xvar <- C(x, rep(0, length(levels(x))), 1) # puts zero in the X matrix 
 attr(xvar, "call") <- substitute( gamlss.random(data[[scall]], z, w))
 attr(xvar, "df") <- df
 attr(xvar, "lambda") <- lambda
 attr(xvar, "start") <- start
 class(xvar) <- c("smooth", class(xvar))
 xvar
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this function is take from Hastie's gam 
# last change MS Thursday, April 25, 2002 at 13:31
gamlss.random <- function(x, y, w) 
{
# local function----------------------------------------------------------------  
    df.inv <- function(n, df, lambda = sum(n)/df - mean(n), iterations = 10)
       { # given df find lambda
        if(df > length(n))
            return(0)
        current.df <- sum(n/(n + lambda))
        if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
            lambda
        else {
            lambda <- exp(log(lambda) + (current.df - df)/(sum((n * lambda)/(n +lambda)^2)))
            Recall(n, df, lambda, iterations - 1)
             }
        }
#-------------------------------------------------------------------------------
regpen <- function(y, w, x, nw, lambda)
{
  beta <- tapply(w * y, x, sum)/(nw + lambda)
    df <- sum(nw[non.zero]/(nw[non.zero] + lambda))
   out <- list(beta=beta, edf=df)     
}
#-------------------------------------------------------------------------------
          df <- as.vector(attr(x,"df"))
      lambda <- as.vector(attr(x, "lambda"))
       start <- as.vector(attr(x, "start"))
          nw <- tapply(w, x, sum) 
    non.zero <- !is.na(nw)
        sig2 <- tau2 <- 0
           N <- sum(w!=0) # DS+FDB 3-2-14   
           n <- length(x)
    # if lamnda is not 
# case 1 if lambda is know    
if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, w, x, nw, lambda)
           fv <- fit$beta[x]        
 } 
if (is.null(df)&&is.null(lambda)) # case 2 if both lambda and df are NULL
 {
   lambda <- start
   for (it in 1:50) 
         {
           fit  <- regpen(y, w, x, nw, lambda) # fit model
         gamma. <- fit$beta  # get the gamma differences
             fv <-  fit$beta[x]            # fitted values
           sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf)# Monday, March 16, 2009 at 20:00 see LNP page 279
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
         # cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
     # assign(startLambdaName, lambda, envir=gamlss.env)
     #cat("lambda",lambda, '\n')
         }
 }
if (!is.null(df)&&is.null(lambda)) # case 3 : if df are required----------------
 {
  lambda <- df.inv(nw[non.zero], df)
   fit <- regpen(y, w, x, nw, lambda)
 }
  #  if(is.null(df))     df <- sum(non.zero)
  #  if(lambda == 0) lambda <- df.inv(nw[non.zero], df)
  #   df <- sum(nw[non.zero]/(nw[non.zero] + lambda))
    beta <- fit$beta #tapply(w * y, x, sum)/(nw + lambda)
      fv <-   fit$beta[x]
     var <- as.vector(w/(nw[x] + lambda))
    residuals <- as.vector(y -fv )
coefSmo <- list(  coef = fit$beta,
                lambda = lambda, 
                   edf = fit$edf, 
                 sigb2 = tau2, 
                 sige2 = sig2,
                sigb = if (is.null(tau2)) NA else sqrt(tau2),
                sige = if (is.null(sig2)) NA else sqrt(sig2),
                    fv = fv,  
                factor = x,
                    se = sqrt(var))
class(coefSmo) <- "random"
    list(x = seq(along = nw), y = fv, residuals = residuals, var = var, 
         nl.df = fit$edf, lambda=lambda, coefSmo=coefSmo) # MS 
}

plot.random <- function(x,...)
{
  plot(levels(x$factor), x$coef, type="h", xlab="levels", ylab="coefficients")
  abline(h=0)
}

coef.random <- function(object, ...)
{
  as.vector(object$coef)
}

fitted.random<- function(object, ...)
{
  as.vector(object$fv)
}

print.random  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("Randon effects fit using the gamlss function random() \n")
#  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
#  cat("Coefficients \n")
#  co <- coef(x)
#  cat("\n  ", names(co), "\n")
#  cc <-simplify2array(co, higher=TRUE)
#  cat(cc, " \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
  cat("Random effect parameter sigma_b:", format(signif(x$sigb)), "\n")  
  cat("Smoothing parameter lambda     :", format(signif(x$lambda)), "\n") 
}
