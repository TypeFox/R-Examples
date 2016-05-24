# The penRegQ function fits a penalsed B-spline models to y a and x
# It is similar to the function PenReg but it uses the Q-function (margional local Likelihood)
# to estimate the sigmas
# created by MS 23--4-12
# TO DO: it needs prediction
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
penRegQ <- function(y, x,
                weights = rep(1,length(y)), # for weighted observations 
                  order = 2,
                  start = 10,
                   plot = FALSE,
                 lambda = NULL,  
                  inter = 20, 
                 degree = 3,
             optim.proc = c("nlminb", "optim"), # ,  "genoud"
          optim.control = NULL)                       
 {
# ---------------------------------------------------
# local functions  
# creates the basis for p-splines
 bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
  {
 tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
# DS xl= min, xr=max, ndx= number of points within 
# Construct B-spline basis
     dx <- (xr - xl) / ndx # DS increment 
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx) 
      P <- outer(x, knots, tpower, deg)# calculate the power in the knots
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
      B <- (-1) ^ (deg + 1) * P %*% t(D) 
      B 
  }
#--------------------------------------------------
#--------------------------------------------------
 # this function from nlme of Pinheiro and Bates 
 HessianPB<-function (pars, fun, ..., .relStep = (.Machine$double.eps)^(1/3), 
                      minAbsPar = 0) 
 {
   pars <- as.numeric(pars)
   npar <- length(pars)
   incr <- ifelse(abs(pars) <= minAbsPar, minAbsPar * .relStep, 
                  abs(pars) * .relStep)
   baseInd <- diag(npar)
   frac <- c(1, incr, incr^2)
   cols <- list(0, baseInd, -baseInd)
   for (i in seq_along(pars)[-npar]) {
     cols <- c(cols, list(baseInd[, i] + baseInd[, -(1:i)]))
     frac <- c(frac, incr[i] * incr[-(1:i)])
   }
   indMat <- do.call("cbind", cols)
   shifted <- pars + incr * indMat
   indMat <- t(indMat)
   Xcols <- list(1, indMat, indMat^2)
   for (i in seq_along(pars)[-npar]) {
     Xcols <- c(Xcols, list(indMat[, i] * indMat[, -(1:i)]))
   }
   coefs <- solve(do.call("cbind", Xcols), apply(shifted, 2, 
                                                 fun, ...))/frac
   Hess <- diag(coefs[1 + npar + seq_along(pars)], ncol = npar)
   Hess[row(Hess) > col(Hess)] <- coefs[-(1:(1 + 2 * npar))]
   list(mean = coefs[1], gradient = coefs[1 + seq_along(pars)], 
        
        Hessian = (Hess + t(Hess)))
 }
 #--------------------------------------------------
 #--------------------------------------------------
# main function starts here
         scall <- deparse(sys.call())
       optim.p <- match.arg(optim.proc)
        natrue <- FALSE 
            if (any(is.na(y)))
           {
                       natrue <- TRUE
                            yy <- y # we need this for plotting
             weights[is.na(y)] <- 0
                   y[is.na(y)] <- 0
           }
       if (any(weights==0))
           {
             natrue <- TRUE
             yy <- ifelse(weights==0, 0, y) 
           } 
    giveLambda <- !is.null(lambda)
        lambda <- if (is.null(lambda)) start else lambda 
             N <- length(y)
         inter <- if (N<100) 10 else inter # this is to prevent singularities when length(x) is small
            xl <- min(x)
            xr <- max(x)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)   
             X <- bbase(x, xmin, xmax, inter, degree) # create the basis
             r <- ncol(X)
             D <- if(order==0) diag(r) else diff(diag(r), diff=order)              
    #        E <- diag.spam(N)
    #        D <- diff(E, diff = order)   
             G <- t(D)%*%D 
            XW <- weights * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + lambda*G, t(XW) %*% y)
            fv <- X %*%beta
            tr <- sum(diag(solve(XWX + G, XWX)))
          sige <- sum(weights*(y-fv)^2)/(N-tr) 
          sigb <- sum((D%*%beta)^2)/(tr)
        params <- c(logsige <- log(sige), logsigb<-log(sigb) )
 if (giveLambda)
{
  fit <- list(fitted = fv, edf = tr, lambda = lambda, order = order, sig2 =  sige, 
                     sigmae = sqrt(sige), tau2 = sigb, sigmab = sqrt(sigb), value.of.Q = NULL,
                        par = params, y = y, weights = weights, N = N,
                      call = scall, rss = sum(weights*(y-fv)^2),
             aic = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE)))+2*(tr) , 
             sbc = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE)))+log(N)*(tr),
             deviance = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE))))
   return(fit)
 }
#-------------------------------------------------------
#-------------------------------------------------------
    Qf<- function(par) 
     { 
  #    cat("par", exp(par), "\n")
     # browser()
       lambda <- exp(par[1]-par[2])
      # if (lambda<=1.0e-7) browser();lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
      # if (lambda==1.0e+7) browser();lambda<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
      #cat("lambda",lambda,"\n")
         beta <- solve(XWX + lambda*G, t(XW) %*% y)
           fv <- X %*%beta  
            b <- D %*% as.vector(beta)
            #browser()
           D3 <- determinant(exp(-par[1])*XWX+exp(-par[2])*G)$modulus
            f <- -(N/2)*log(2*pi*exp(par[1]))+              #1
                 .5*sum(log(weights))-   #2 this do not depends on parameters
                 sum(weights*(y-fv)^2)/(2*exp(par[1]))+     #3
                  (r/2)*log(2*pi)-                          #4 
                 ((r-order)/2)*log(2*pi*exp(par[2]))-       #5
                 sum(b^2)/(2*exp(par[2]))-                  #6
                   .5*D3                                    #7
        #   cat("f", -f, "\n")      
           -f
     }
 #--------------------------------------------------------
 #--------------------------------------------------------
 switch(optim.p, "nlminb" = { # this do not allow se's'
            out <- nlminb(start = params, objective = Qf, lower = c(-16, -16),  upper = c(16, 16), control=optim.control)
            if (out$convergence > 0) # I took this from Ripley
        warning("possible convergence problem: optim gave code=", 
                out$convergence, " ", out$message)
    out$hessian <- optimHess(out$par, Qf)  
    #out$hessian <- HessianPB(pars=out$par, fun=Qf)$Hessian
     value.of.Q <- out$objective
                         },
                 "optim"={
           out <- optim(par= params, fn = Qf, lower =  c(-16, -16), upper = c(16, 16),method= "L-BFGS-B", control=optim.control)
      if (out$convergence > 0) # I took this from Ripley
        warning("possible convergence problem: optim gave code=", 
                out$convergence, " ", out$message)
    out$hessian <- optimHess(out$par, Qf)              
    #out$hessian <- HessianPB(pars=out$par, fun=Qf)$Hessian
                  # geting the hessian using nlme    fdHess(out$par, Qf)
            value.of.Q <- out$value
                           },
#                  "genoud"={
#                    B <- matrix(c(-16,-16, 16, 16), ncol=2)
#                    # spead is INCREASED  by reducing  pop.size and increasing BFGSburnin
#              out <- genoud( Qf, nvars=2, max=FALSE, pop.size=30,  Domains=B, boundary.enforcement=2, #gradient.check=TRUE,BFGSburnin=10, max.generations=30, print.level=0)#  hessian=TRUE
#      out$hessian <- optimHess(out$par, Qf)  
#     # out$hessian <- HessianPB(pars=out$par, fun=Qf)$Hessian     
#            # out <- genoud( Qf, nvars=2, max=FALSE, pop.size=10,  Domains=B, boundary.enforcement=2, #hessian=TRUE, BFGSburnin=50, P1=0, P2=50, P3=50, P4=0, P5=0, P6=0, P7=0, P8=0, P9=0)
#            #out <- genoud( Qf, nvars=2, max=FALSE, pop.size=5000, Domains=B, boundary.enforcement=2, #hessian=TRUE)
#                    value.of.Q <- out$value   
#                    }
    )    
  #   out1<- optim(par= params,            fn = Qf, lower = c(1e-10, 1e-10),  upper = c(Inf, Inf),   method= "L-BFGS-B")
 #  B <- matrix(c(1e-10,1e-10, Inf, Inf), ncol=2)  
 # B <- matrix(c(0.001,0.001, 10000, 10000), ncol=2)
 #    genoud( Qf, nvars=2, max=FALSE, pop.size=3000, Domains=B)
      sige <- exp(out$par[1])
      sigb <- exp(out$par[2])
      shes<- try(solve(out$hessian))
      se <- if (any(class(shes)%in%"try-error")) rep(NA_real_,2) else sqrt(diag(shes))
names(out$par) <- c("log(sige^2)", "log(sigb^2)")

    lambda <- sige/sigb
      beta <- solve(XWX + lambda*G, t(XW) %*% y)
        fv <- X %*%beta  
          b <- D %*% as.vector(beta)     
       tr1 <- order + sum(b^2)/(sigb)
       tr2 <- N-(sum(weights*(y-fv)^2))/(sige)
       
    if (plot) {if (natrue) plot(x, yy, col = gray(0.7))
             else plot(x, y, col = gray(0.7))
               lines(fv[order(x)]~x[order(x)], col = 'red', lwd = 2)}
    # get the output 
           fit <- list(
                fitted.values = as.vector(fv),
                 coefficients = as.vector(beta),
                          edf = tr1,
                          edf2= tr2,
                       lambda = lambda, 
                        order = order, 
                         sig2 =  sige,
                        sigma = sqrt(sige),
                   value.of.Q = value.of.Q,
                         tau2 = sigb,
                       sigmab = sqrt(sigb),
                   value.of.Q = value.of.Q,
                          par = list(par=out$par,se=se),
                         vcov = shes, 
                            y = y,
                      weights = weights,
                            N = N,
                         call = scall,
                       method = "Q function",
                          rss = sum(weights*(y-fv)^2),
                          aic = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE)))+2*(tr1+1), 
                          sbc = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE)))+log(N)*(tr1+1),
                     deviance = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE))))
 class(fit) <- c("penRegQ", "penReg")
  fit
  }
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# methods 
coef.penRegQ<-function(object, type = c("sigmas", "betas"), ...) 
{
  type <- match.arg(type)
  if (type == "betas") return(object$coefficients)
  else 
  {
    coef <- object$par[["par"]]
    attr(coef, "se") <-  object$par[["se"]]
    return(coef) 
  }
  coef <- object$par[["par"]]
}
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
deviance.penRegQ<-function(object,type = c("GD", "Marginal"), ...) 
{
  type <- match.arg(type)
  if (type == "GD") return(object$deviance)
  else return(2*object$value.of.Q)
}
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
print.penRegQ  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("\nPenalised Regression fit")
  cat("Fitting method:", deparse(x$method), "\n")
  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
  cat("Coefficients")
  co <- coef(x)
  cat("\n  ", names(co), "\n")
  cc <-simplify2array(co, higher=TRUE)
  cat(cc, " \n")
  cat("\n Degrees of Freedom for the fit:", x$df, "Residual Deg. of Freedom  ", 
      x$N-x$df, "\n")
  cat("Global Deviance:    ", format(signif(x$deviance)), 
      "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
      format(signif(x$sbc)), "\n")
  invisible(x)
}
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
summary.penRegQ  <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("\nPenalised regression fit","\n")
  cat("Fitting method:", deparse(object$method), "\n")
  cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
  coef  <- coef(object)
  se.coef <- attributes(coef)$se
  tval <- coef/se.coef
  matcoef <- cbind(coef, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) <- list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  
  cat("\n Degrees of Freedom for the fit:", object$df, "Residual Deg. of Freedom  ", 
      object$N-object$df, "\n")
  cat("Global Deviance:    ", format(signif(object$deviance)), 
      "\n            AIC:    ", format(signif(object$aic)), "\n            SBC:    ", 
      format(signif(object$sbc)), "\n")
  invisible(object)
}
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

