################################################################################
#### AUTHOR:    Arnost Komarek                                              ####
####            29/04/2004                                                  ####
####                                                                        ####
#### FILE:      smoothSurvReg.R                                             ####
####                                                                        ####
#### FUNCTIONS: smoothSurvReg                                               ####
####            dextreme                                                    ####
####            dstextreme                                                  ####
####            dstlogis                                                    ####
####            piece                                                       ####
####                                                                        ####
#### CHANGES:  20140424 small piece of code added to treat non-zero fail    ####
####                    when doing the grid search for optimal lambda       ####
####                                                                        ####
################################################################################

### ====================================================================================
### smoothSurvReg: Survival regression with smoothed error distribution (main function)
### ====================================================================================
## formula
## data
## subset
## na.action ... na.fail is default, it is not recommended to change it when logscale
##               depends on covariates
## init.beta ... c(initial intercept, initial betas)
##                  give NA's for values whose initials are to be found automatically
## init.scale ... initial value for scale
## init.c ....... initial values for c coefficients
##                  vector of length nsplines with all components 0 < c_i < 1
##                  which sums up to 1
## init.dist .... preferable distribution used in 'survreg' to find initial values
##                if "best" the best fitting distribution is used
##                rayleigh and exponential are changed into weibull
## aic .......... should I search for the "best" lambda using AIC?
## lambda .. grid of lambdas to be searched for the best AIC
##                (I recommend to start with bigger lambdas)
##                if aic = FALSE, only the first lambda is used
## model
## control
## ...
smoothSurvReg <- function(formula = formula(data),
                          logscale = ~1,
                          data = parent.frame(),
                          subset,
                          na.action = na.fail,
                          init.beta,
                          init.logscale,
                          init.c,
                          init.dist = "best",
                          update.init = TRUE,
                          aic = TRUE,
                          lambda = exp(2:(-9)),
                          model = FALSE,
                          control = smoothSurvReg.control(),
                          ...)
{
   ## Give a list of control values
   ## for the fitting process.
   if (missing(control)) control <- smoothSurvReg.control(...)

   ## Load survival package if not loaded
   #mamho <- require(survival)
   #if (!mamho)
   #   stop("I need 'survival' package to be installed. ")

   ## Check initial distribution
   dist.allowed <- c("lognormal", "loggaussian", "loglogistic", "weibull", "rayleigh", "exponential", "best")
   dist.now <- pmatch(init.dist, dist.allowed, nomatch = 0)
   if (dist.now == 0){
      stop("Unknown initial distribution. ")
   }

   ## Give a function call to be recorded in a resulting object.
   call <- match.call(expand.dots = TRUE)

   ## Give a function call to be work with.
   m <- match.call(expand.dots = FALSE)

   ## 'survreg' to compute reference values of the loglikelihood
   ## and to find initial estimates
   ## It will also check many things for consistency
   ##    loglik.ref = max(loglikelihoods of the three fitted models)
   temp <- c("", "formula", "data", "subset", "na.action")
   fit.logn <- m[match(temp, names(m), nomatch=0)]
   fit.logl <- m[match(temp, names(m), nomatch=0)]
   fit.weib <- m[match(temp, names(m), nomatch=0)]
   fit.logn[[1]] <- as.name("survreg")
   fit.logl[[1]] <- as.name("survreg")
   fit.weib[[1]] <- as.name("survreg")
   fit.logn$dist <- "lognormal"          
   fit.logl$dist <- "loglogistic"
   fit.weib$dist <- "weibull"
   #fit.logn$failure <- 2       ### this argument used to be passed to survreg.control in R 2.8.1 and older
   #fit.logl$failure <- 2       ### * starting from R 2.9.0, 'failure' is no more argument of survreg.control
   #fit.weib$failure <- 2
   fit.logn <- eval(fit.logn, parent.frame())
   fit.logl <- eval(fit.logl, parent.frame())
   fit.weib <- eval(fit.weib, parent.frame())
   loglik.three <- numeric()
   fit0 <- NULL

   if (is.null(fit.logn$fail)){
      loglik.three <- c(loglik.three, fit.logn$loglik[2])
      fit0 <- fit.logn
   }
   else
      loglik.three <- c(loglik.three, -Inf)

   if (is.null(fit.logl$fail)){
      loglik.three <- c(loglik.three, fit.logl$loglik[2])
      fit0 <- fit.logl
   }
   else
      loglik.three <- c(loglik.three, -Inf)

   if (is.null(fit.weib$fail)){
      loglik.three <- c(loglik.three, fit.weib$loglik[2])
      fit0 <- fit.weib
   }
   else
      loglik.three <- c(loglik.three, -Inf)

   if (is.null(fit0))     ## none of the three 'survreg' distributions was successful
      stop("Sorry but neither 'survreg' is able to fit the model. ")

   dist.user <- switch(dist.now, 1, 1, 2, 3, 3, 3, 4)   ## 1 = lognormal, 2 = loglogistic, 3 = weibull, 4 = best
   loglik.ref <- max(loglik.three)                      ## It must be finite value now (due to the previous rows)
   dist.best <- which.max(loglik.three)
   if (dist.user < 4){
       loglik.user <- loglik.three[dist.user]
       if (loglik.user == -Inf){
           dist.user <- dist.best
           loglik.user <- loglik.ref
       }
   }
   else{
       dist.user <- dist.best
       loglik.user <- loglik.ref
   }

   fit0 <- switch(dist.user, fit.logn, fit.logl, fit.weib)
   init.dist <- switch(dist.user, "lognormal", "loglogistic", "weibull")

      ### for compatibility with the folowing code which is older
   dist.now <- switch(dist.user, 1, 3, 4)       ## 1 = lognormal, 3 = loglogistic, 4 = weibull

   ## Which of the following formal argumets were really used in a
   ## function call?
   ## Store in m only these, throw away remaining ones.
   ## "" states actually for a name of the function.
   m.keep <- m
   temp <- c("", "formula", "data", "subset", "na.action")
   m <- m[match(temp, names(m), nomatch=0)]

   ## Change the value of m[[1]] from "survreg" into "model.frame".
   m[[1]] <- as.name("model.frame")

   ## Which functions should be considered to be special
   ## when constructing a terms object from a formula.
   special <- c("strata", "cluster", "frailty")

   ## Construct a terms object from a formula.
   Terms <- if(missing(data)) terms(formula, special)
            else              terms(formula, special, data=data)

   ## Neither strata nor cluster nor frailties are allowed.
   if(!is.null(attr(Terms,"specials")$strata)){
      stop("Strata in a model formula not implemented for this function. ")
   }
   if(!is.null(attr(Terms,"specials")$cluster)){
       stop("Cluster in a model formula not implemented for this function. ")
   }
   if(!is.null(attr(Terms,"specials")$frailty)){
       stop("Frailty in a model formula not implemented for this function. ")
   }

   is.intercept <- ifelse(attr(Terms,"intercept") == 1, TRUE, FALSE)

   ## Change the formula part of m object into
   ## somewhat more complex object with class terms.
   m$formula <- Terms

   ## Evaluate m. At this moment, m is something like
   ## model.frame(formula=Surv(x,event)~cov1, data=data etc.).
   ## m has now mode "call".
   m <- eval(m, parent.frame())
        ### The mode of m is now "list". But it has also
        ### many useful attributes containing lots of information.

   ## Extract the response.
   ## (it is still survival object)
   Y <- model.extract(m, "response")
   if (!inherits(Y, "Surv"))
      stop("Response must be a survival object. ")

   ## Create a design matrix.
   X <- model.matrix(Terms, m)
   n <- nrow(X)
   nvar <- ncol(X)
   if (nvar <= 0)
      stop("Invalid design matrix. ")

   ## Check whether the response does not have type 'counting'
   ## which is not allowed by smoothSurvReg function.
   type <- attr(Y, "type")
   if (type== 'counting') stop ("Invalid survival type ('counting' is not implemented). ")

   ## Create an offset term (if not presented give all 0 into it).
   offset <- attr(Terms, "offset")
   if (!is.null(offset)) offset <- as.numeric(m[[offset]])
   else                  offset <- rep(0, n)

   ## Log-transformation of the response
     tranfun <- function(y) log(y)
     dtrans <- function(y) 1/y
     exactsurv <- (Y[,ncol(Y)] == 1)   ## rows with exact survival

     ## For exact survivals, log(jacobian) has to be added to the loglikelihood.
     ## (since it will be further worked with transformed variable)
     if (any(exactsurv)) logcorrect <- sum(log(dtrans(Y[exactsurv,1])))
     else                logcorrect <- 0

     ## Transform it
     if (type == 'interval') {
        if (any(Y[,3]==3)) Y <- cbind(tranfun(Y[,1:2]), Y[,3])
        else               Y <- cbind(tranfun(Y[,1]), Y[,3])
     }
     else if (type=='left'){
             Y <- cbind(tranfun(Y[,1]), 2-Y[,2])   ## change 0 indicator into 2 indicating left censoring
          }
          else  ## type = 'right' or 'interval2'
             Y <- cbind(tranfun(Y[,1]), Y[,2])

     if (!all(is.finite(Y))) stop("Invalid survival times for this distribution (infinity on log-scale not allowed). ")

   
   ## Design matrix for logscale
     common.logscale <- FALSE                            ## indicator whether only intercept is included in log(scale) formula
     if (!control$est.scale) common.logscale <- TRUE
     else{
       if (!match("logscale", names(m.keep), nomatch = 0)) common.logscale <- TRUE
       else{
         tempR <- c("", "logscale", "data", "subset", "na.action")
         mR <- m.keep[match(tempR, names(m.keep), nomatch=0)]
         mR[[1]] <- as.name("model.frame")
         names(mR)[2] <- "formula"
         TermsR <- if(missing(data)) terms(logscale)
                   else              terms(logscale, data = data)
         lTR <- length(attr(TermsR, "variables"))
         if (lTR == 1 & !attr(TermsR, "intercept")){        ## nothing specified, include at least intercept
           attr(TermsR, "intercept") <- 1
           common.logscale <- TRUE
         }
         else{
           if (lTR == 1 & attr(TermsR, "intercept")){        ## the only term is the intercept
             common.logscale <- TRUE
           }
           else{
             mR$formula <- TermsR
             mR <- eval(mR, parent.frame())
             if (attr(TermsR, "intercept")){    ## the intercept in included
               ## HERE: DO NOTHING
             }
             Z <- model.matrix(TermsR, mR)
           }
         }
       }
     }       
     if (common.logscale){
       Z <- matrix(rep(1, n), ncol = 1)
       colnames(Z) <- "(Intercept)"
     }
     names.logscale <- colnames(Z)   
     n.logscale <- length(names.logscale)     

   
## Initial values for BETA coefficients and the INTERCEPT
## ------------------------------------------------------
   ninit.beta <- dim(X)[2]      ## number of initial values for beta parameters

      # All initial values from survreg
   if (missing(init.beta) || is.null(init.beta)){
       init.beta <- fit0$coefficients
       if (is.intercept){
          if (dist.now %in% 1:3)      ## lognormal or loglogistic initial distribution
            init.beta[1] <- init.beta[1]
          else
            if (dist.now %in% 4:6)    ## weibull initial distribution
                 init.beta[1] <- init.beta[1] - 0.5772*fit0$scale
            else
                 stop("Unknown initial distribution. ")
       }
   }

      # (Some) initial values from the user
   else{
       if (length(init.beta) != ninit.beta) stop("Invalid length of the vector 'init.beta'. ")

         # Intercept
       if (is.na(init.beta[1]) && is.intercept)
          if (dist.now %in% 1:3){      ## lognormal or loglogistic initial distribution
            init.beta[1] <- fit0$coefficients[1]
          }
          else
            if (dist.now %in% 4:6){ ## weibull initial distribution
                 init.beta[1] <- fit0$coefficients[1] - 0.5772*fit0$scale
            }
            else
                 stop("Unknown initial distribution. ")

         # Betas (except the intercept)
       first.nonintercept <- ifelse(is.intercept, 2, 1)
       ninit.betareal <- ifelse(is.intercept, ninit.beta-1, ninit.beta)
       if ((ninit.beta > 1 && is.intercept) || (ninit.beta == 1 && !is.intercept)){
          betainit <- init.beta[first.nonintercept:ninit.beta]
          betafit0 <- fit0$coefficients[first.nonintercept:ninit.beta]
          betainit[is.na(betainit)] <- betafit0[is.na(betainit)]
          init.beta[first.nonintercept:ninit.beta] <- betainit
       }

       if (is.null(names(init.beta))){
           namesx <- dimnames(X)[[2]]
           if (is.intercept)
               if (is.null(namesx))
                   namesx <- c("(Intercept)", paste("beta", 1:ninit.betareal, sep=""))
               else
                   namesx[1] <- "(Intercept)"
           else
               if (is.null(namesx))
                   namesx <- paste("beta", 1:ninit.betareal, sep="")
           names(init.beta) <- namesx
       }
   }

### Initial values for log(SCALE) parameters
### ----------------------------------------
      # Initial value from survreg
   if (fit0$scale <= 0)
      temp.logscale <- log(0.01)
   else
      if (dist.now %in% 1:2)           ## lognormal initial distribution
        temp.logscale <- log(fit0$scale)
      else
        if (dist.now == 3)             ## loglogistic initial distribution
           temp.logscale <- log(fit0$scale * (pi/sqrt(3)))
        else
           if (dist.now %in% 4:6)     ## weibull initial distribution
              temp.logscale <- log(fit0$scale * (pi/sqrt(6)))
           else
              stop("Unknown initial distribution ")
   
   if (missing(init.logscale) || is.null(init.logscale) || sum(is.na(init.logscale))){       
       init.logscale <- c(temp.logscale, rep(0, n.logscale - 1))
   }   

      # Initial value from the user
   else{
       if (length(init.logscale) != n.logscale) init.logscale <- c(temp.logscale, rep(0, n.logscale - 1))
       else                                     init.logscale <- init.logscale
   }
   names(init.logscale) <- names.logscale

## Initial values for G-SPLINE coefficients
## (if given by the user and est.c I use only the first g-3 ones
##  the rest is calculated from these at the beginning)
## --------------------------------------------------------------
   if (control$est.c){     ## nsplines is also at least 4

      if (missing(init.c) || is.null(init.c)){
             ## try to approximate the "best" distribution according to 'survreg'
             ## or the distribution required by the user
             best.dens <- switch(dist.user, "dnorm", "dstlogis", "dstextreme")
             init.c <- find.c(control$knots, control$sdspline, best.dens)
             if (!is.null(init.c)){
                i1 <- which.max(init.c)
                i2 <- which.max(init.c[-i1]); i2 <- ifelse(i2 < i1, i2, i2 + 1)
                i3 <- which.max(init.c[-c(i1, i2)]); i3 <- ifelse(i3 < min(i1, i2), i3, ifelse(i3 < max(i1, i2) - 1, i3 + 1, i3 + 2))
                last.three.temp <- c(i1, i2, i3)
                init.c <- give.c(control$knots, control$sdspline, last.three.temp, init.c[-last.three.temp])
                init.c[init.c < 1e-5] <- 1e-5
             }
             else{
                ## USE ANOTHER METHOD ==> LATER ON (MAYBE)
                stop("Singularity when computing initial c's, try to give your own initial c's or a's  ")
             }
      }
      else{
              if(length(init.c) != control$nsplines) stop("Incorrect length of the vector of initial c coefficients. ")
              if((sum(init.c) < 0.99) || (sum(init.c) > 1.01)) stop("Sum of initial c coefficients is not 1. ")
              if((sum(init.c < 0) > 0) || (sum(init.c > 1) > 0)) stop("All c coefficients must be between 0 and 1. ")
              init.c <- give.c(control$knots, control$sdspline, control$last.three, init.c[-control$last.three])
              init.c[init.c < 1e-5] <- 1e-5
       }
   }

   else{          ## c's are not estimated
      if (missing(init.c) || is.null(init.c))
            if (control$nsplines == 1) init.c <- 1
            else                       stop("Initial a's or c's must be given. ")
      else{
            if(length(init.c) != control$nsplines) stop("Incorrect length of the vector of initial c coefficients. ")
            if((sum(init.c) < 0.99) || (sum(init.c) > 1.01)) stop("Sum of initial c coefficients is not 1. ")
            if((sum(init.c <= 0) > 0) || (sum(init.c > 1) > 0)) stop("All c coefficients must be between 0 and 1. ")
       }
   }

## Put all initials into a list
   initials <- list(beta = init.beta, logscale = init.logscale, ccoef = init.c)

   if (control$debug == 1){
      cat("\ncontrol:\n"); print(control)
      cat("\ninitials:\n"); print(initials)
      cat("\n"); print(summary(fit0))
   }

   
## Do not use searching for the best lambda if !est.c or if maxiter == 0
## ---------------------------------------------------------------------
   if (!control$est.c)         aic <- FALSE     ## all real lambda's give same fit
   if (control$maxiter == 0)   aic <- FALSE     ## user wants the derivatives at some point
   if (!aic)                   lambda <- lambda[1]
   initials.aic <- initials

## Search for the best AIC if desired (otherwise fit it only once for the first lambda)
   lambda <- lambda[order(lambda, decreasing = TRUE)]
   nlam <- length(lambda)
   if (sum(lambda < 0) > 0) stop("'lambda' must contain only non-negative values. ")
   fit.aic <- list()
   aic.values <- numeric()
   df.previous <- -1e40
   fail.previous <- 0
   df.values <- numeric()
   pll.values <- numeric()
   ll.values <- numeric()
#   df2.values <- numeric()
   nofpar.values <- numeric()
   problem.look <- numeric()
   problem <- numeric()
   search <- TRUE
   m <- 1
   warn.opt <- options("warn")$warn
   options(warn = -1)
   while (search){
      control$lambda.use <- lambda[m]
      if (control$info){
         cat("\n\n=================================================")
      }
      cat("\nFit with Log(Lambda) = ", log(lambda[m]), sep="")
      if (!control$info) cat(",  ")
      
      fitA <- smoothSurvReg.fit(X, Z, Y, offset, correctlik = logcorrect,
                                init = initials.aic, controlvals = control, common.logscale = common.logscale)

      ### === Code added on 20140424 === ###
      if (fitA$fail & m > 1){       ### Try also original initials if any type of fail
        fitAB <- smoothSurvReg.fit(X, Z, Y, offset, correctlik = logcorrect,
                                   init = initials, controlvals = control, common.logscale = common.logscale)
        if (fitAB$fail == 0) fitA <- fitAB
      }  
      ### === End of code added on 20140424 === ###
        
      fit.aic[[m]] <- fitA

      if (fitA$fail >= 99){
         fitA$aic <- NA
         fitA$degree.smooth$df <- NA
         fitA$loglik["Penalized Log Likelihood"] <- NA
         fitA$loglik["Log Likelihood"] <- NA
         fitA$degree.smooth[["Number of parameters"]] <- NA
         fitA$iter <- NA
      }
      else{
         sd.regres <- as.numeric(fitA$regres[["Std.Error"]])
         sd2.regres <- as.numeric(fitA$regres[["Std.Error2"]])
         sd.nans <- sum(is.na(sd.regres))
         sd2.nans <- sum(is.na(sd2.regres))
         if (n.logscale <= 1){
           if(control$est.scale && sd.nans >= 2) fitA$fail <- fitA$fail + 40
           if(!control$est.scale && sd.nans >= 1) fitA$fail <- fitA$fail + 40
         }
         else{
           if(control$est.scale && sd.nans >= 1) fitA$fail <- fitA$fail + 40
           if(!control$est.scale && sd.nans >= 0) fitA$fail <- fitA$fail + 40
         }           
      }

      aic.values <- c(aic.values, fitA$aic)
      df.values <- c(df.values, fitA$degree.smooth$df)
#      df2.values <- c(df2.values, fitA$degree.smooth$df2)
      pll.values <- c(pll.values, as.numeric(fitA$loglik["Penalized Log Likelihood"]))
      ll.values <- c(ll.values, as.numeric(fitA$loglik["Log Likelihood"]))
      nofpar.values <- c(nofpar.values, fitA$degree.smooth[["Number of parameters"]])
      problem <- c(problem, fitA$fail)
      problem.look <- c(problem.look, fitA$fail)

      cat("AIC(", lambda[m], ") = ", fitA$aic, sep="")
      cat(",  df(", lambda[m], ") = ", fitA$degree.smooth$df, sep="")
#      cat(",  df2(", lambda[m], ") = ", fitA$degree.smooth$df2, sep="")
#      cat(",  n of param.(", lambda[m], ") = ", fitA$degree.smooth[["Number of parameters"]], sep="")
      cat(",  ", fitA$iter, " iterations", sep = "")
      cat(",  fail = ", fitA$fail, sep = "")
      if (m == length(lambda)) cat("\n")

      if (fitA$fail < 99){
         if (update.init && fitA$fail == 0 && fitA$degree.smooth$df > df.previous){   ## update the initials
             initials.aic$beta <- as.numeric(fitA$regres[1:ninit.beta, 1])
             if (control$est.scale){
               if (n.logscale == 1) initials.aic$logscale <- as.numeric(fitA$regres[ninit.beta+2, 1])
               else                 initials.aic$logscale <- as.numeric(fitA$regres[(ninit.beta+1):(ninit.beta+n.logscale), 1])
             }               
             if (control$est.c)     initials.aic$ccoef <- as.numeric(fitA$spline[, "c coef."])
         }

         if (fitA$fail == 0 && fitA$degree.smooth$df < df.previous - 1 && fail.previous == 0)
             problem.look[length(problem.look)] <- 99
         else
             df.previous <- fitA$degree.smooth$df
      }
      fail.previous <- fitA$fail

      m <- m + 1
      if (m > nlam) search <- FALSE
   }    ## end of 'while (search)'
   options(warn = warn.opt)

## If all lambdas caused problems then the last fit with fail < 99 is presented
##   otherwise the fit with fail == 0 and highest AIC is presented as the last one
   n.aic <- length(aic.values)
   if (sum(problem.look >= 99) == n.aic){
      fit.smooth <- fit.aic[[1]]    ## all fits are same so that I give that first one
   }
   else{
      if (sum(problem.look >= 20) == n.aic){  ## all fits are without df
          give <- max((1:n.aic)[problem.look < 99])   ## return the last one with fail < 99
      }
      else{
         if (sum(problem.look > 0) == n.aic){  ## all fits are problematic but at least one of them has fail < 99 and fail < 20
             give <- max((1:n.aic)[problem.look < 20])   ## return the last one with fail < 99 and fail < 30
         }
         else{   ## there is at least one fit with fail == 0
             give <- (1:n.aic)[problem.look == 0]
             mat.help <- rbind(give, aic.values[problem.look == 0])
             give <- which.max(aic.values[problem.look == 0])
             give <- mat.help[1, give]
         }
      }
      fit.smooth <- fit.aic[[give]]
      control$lambda.use <- lambda[give]
   }

   aic.values <- c(aic.values)
   df.values <- c(df.values)
   nofpar.values <- c(nofpar.values)
   searched <- data.frame(Lambda = lambda, LogLambda = log(lambda), AIC = aic.values, df = df.values,
                     PenalLogLik = pll.values, LogLik = ll.values,
                     nOfParm = nofpar.values, fail = problem)
   colnames(searched)[2] <- "Log(Lambda)"
#   searched <- data.frame(Lambda = lambda, AIC = aic.values, df = df.values, df2 = df2.values,
#                     nOfParm = nofpar.values, fail = problem)


## No fit produced
   if (fit.smooth$fail >= 99){
      warning("No fit is produced ")
      class(fit.smooth) <- 'smoothSurvReg'
      return(fit.smooth)
   }

## Handle possible warnings
   if (fit.smooth$fail > 0){
       for (i in 1:3)
          if (fit.smooth$warning[i, 1] != "OK") warning(fit.smooth$warning[i, 1])
   }

## Further manipulation with the results
   na.action <- attr(m, "na.action")
   if (length(na.action)) fit.smooth$na.action <- na.action
   fit.smooth$terms <- Terms
   fit.smooth$formula <- as.vector(attr(Terms, "formula"))
   fit.smooth$call <- call
   fit.smooth$init.dist <- init.dist
   if (model) fit.smooth$model <- m
   fit.smooth$x <- X
   fit.smooth$y <- Y
   fit.smooth$z <- Z

   ## Initial c coefficients
   knotname <- paste("knot[",1:control$nsplines,"]", sep = "")
   sd.spline <- rep(control$sdspline, control$nsplines)
   fit.smooth$init.spline <- data.frame(Knot = as.numeric(control$knots),
                                        SD.spline = as.numeric(sd.spline),
                                        c.coef = as.numeric(initials$ccoef)
                             )
   rownames(fit.smooth$init.spline) <- knotname
   colnames(fit.smooth$init.spline) <- c("Knot", "SD basis", "c coef.")

   
   ## Put initial alpha, beta and log(scale) pars. estimates into the resulting object
   if (common.logscale){
     fit.smooth$init.regres <- c(initials$beta, initials$logscale, exp(initials$logscale))
     names(fit.smooth$init.regres) <- c(names(initials$beta), "Log(scale)", "Scale")
   }
   else{
     fit.smooth$init.regres <- c(initials$beta, initials$logscale)
     names(fit.smooth$init.regres) <- c(names(initials$beta), paste("LScale.", names.logscale, sep = ""))     
   }     
   fit.smooth$init.regres <- data.frame(Value = fit.smooth$init.regres)

   ## Compute mean and variance of the error distribution
   ## This should be zero and one if c's are estimated
   ccoef <- fit.smooth$spline[["c coef."]]
   mean.error <- sum(ccoef * control$knots)
   var.error <- sum(ccoef * (sd.spline^2 + (control$knots)^2)) - mean.error^2
   sd.error <- ifelse(var.error >=0 , sqrt(var.error), NaN)

   ## Compute adjusted intercept and scale
   ## (after taking into account mean and scale of the error term)
   if (is.intercept)
      mu0 <- fit.smooth$regres["(Intercept)","Value"]
   else
      mu0 <- 0

   if (control$est.scale)
     if(common.logscale) s0 <- fit.smooth$regres["Scale","Value"]
     else                s0 <- NA
   else
     s0 <- fit.smooth$init.regres["Scale", "Value"]

   if (common.logscale){
     int.adj <- mu0 + (s0 * mean.error)
     scale.adj <- s0 * sd.error
   }
   else{
     int.adj <- NA
     scale.adj <- NA
   }

   fit.smooth$adjust <- data.frame(Value = c(int.adj, scale.adj))
   rownames(fit.smooth$adjust) <- c("(Intercept)", "Scale")

   fit.smooth$error.dist <- data.frame(Mean = mean.error, Var = var.error, SD = sd.error)
   rownames(fit.smooth$error.dist) <- "Error distribution:  "

   fit.smooth$searched <- searched

   ## Add indicator of the presence of the intercept in the model to estimated
   fit.smooth$estimated <- c(is.intercept, fit.smooth$estimated)
   names(fit.smooth$estimated) <- c("(Intercept)", "Scale", "ccoef", "common.logscale")

   class(fit.smooth) <- 'smoothSurvReg'
   return(fit.smooth)
}


### =============================================
### piece: Evaluate a piecewise constant function
### =============================================
## (used when computing initial c coefficients)
## Assumption: breaks are sorted
piece <- function(x, breaks, values){
   x <- x[order(x)]
   if (length(breaks) != (length(values)+1))
      stop("Badly defined piecewise constant function ")
   fx <- numeric(length(x))

   fx[x<=breaks[1]] <- 0
   for (i in 1:(length(breaks)-1))
      fx[breaks[i] < x & x <= breaks[i+1]] <- values[i]
   fx[x > breaks[length(breaks)]] <- 0

   return(fx)
}


### ================================================
### dextreme: Density of extreme value distribution
### ================================================
dextreme <- function(x, alpha=0, beta=1){
  value <- (1/beta)*exp((x-alpha)/beta)*exp(-exp((x-alpha)/beta))
  return(value)
}


### ============================================================================
### dstextreme: Density of standardized (E=0, var=1) extreme value distribution
### ============================================================================
dstextreme <- function(x){
  beta <- sqrt(6)/pi
  alpha <- beta*0.5772
  value <- dextreme(x, alpha, beta)
  return(value)
}


### ========================================================
### dstlogis: Density of standardized logistic distribution
### ========================================================
dstlogis <- function(x){
  scale <- sqrt(3)/pi
  value <- dlogis(x, 0, scale)
  return(value)
}

