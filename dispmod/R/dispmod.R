##
## Normal dispersion model
##

lm.disp <- function(formula, var.formula=NULL, data = list(), maxit = 30, epsilon = glm.control()$epsilon, subset, na.action = na.omit, contrasts = NULL, offset = NULL)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$method <- mf$contrasts <- mf$... <- NULL
  mf[[1]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE
  mf$epsilon <- mf$maxit <- mf$weights <- NULL
  
  # Defines a model frame for the mean and one for the variance
  if (!is.null(var.formula))
     {  mf.var <- mf[-(which(names(mf)=="formula"))]
            names(mf.var)[2] <- "formula"
                mf     <- mf[-(which(names(mf)=="var.formula"))] }
  else
     { mf.var <- mf }

  # remove response term in variance formula if any
  f <- paste(mf.var$"formula")
  if (length(f)==3)
         mf.var$"formula" <- as.formula(f[-2])
         
  mean.formula <- mf$"formula"   
  var.formula  <- mf.var$"formula"       
  
  ## collect information for the mean function model
  ##
  mf <- eval(mf, sys.frame(sys.parent()))
  mt <- attr(mf, "terms")
  xvars <- as.character(attr(mt, "variables"))[-1]
  if ((yvar <- attr(mt, "response")) > 0) 
      xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) 
             { xlev <- lapply(mf[xvars], levels)
               xlev[!sapply(xlev, is.null)] }
  y <- model.response(mf, "numeric")
  n <- NROW(y)
  offset <- model.offset(mf)
  if (!is.null(offset) && length(offset) != n) 
      stop(paste("Length of the offset", length(offset), ", must equal", 
                 n, " the number of cases"))
  if (is.empty.model(mt)) 
      stop(paste("No model specified!"))
  # mean function predictors
  x <- model.matrix(mt, mf, contrasts)

  ## collect information for the variance function model
  ##
  mf.var <- eval(mf.var, sys.frame(sys.parent()))
  mt.var <- attr(mf.var, "terms")
  if (!is.null(model.extract(mf.var, "response")))
     warning("response ignored in the formula for the variance function!")
  # variance function predictors
  z <- model.matrix(mt.var, mf.var, contrasts)

  # function to the deviance for dispersion models.
  disp.dev <- function(mod, mod.var) 
  { 
        r <- mod$residuals
        dev <- mod$deviance
        n <- length(r)
        s <- dev/n * mod.var$fitted.values
        sum( log(2*pi) + log(s) + r^2/s)
  }

  ##  Start iterative procedure  ##
  
  # estimate mean function via OLS.
  mod     <- glm.fit(x, y, family=gaussian(identity))
  # estimate variance function 
  mod.var <- glm.fit(z, mod$residuals^2, family=Gamma(log))

  initial.dev <- n * (log(2*pi) + 1 + log(sum(mod$residuals^2)/n))
  dev0 <- Inf
  dev1 <- initial.dev

  # loop until change in deviance is less negligible
  i <- 1
  cat("Iteration ", i, ": deviance = ", format(dev1), "\n", sep="")
  while( abs(dev1-dev0) > epsilon )
  {
        if (i > maxit) 
       { warning("algoritm not converged after ", i, " iterations!")
         return() }

        i <- i + 1
        # working weights
        w <- 1/(mod.var$fitted.values)
        # re-estimate mean function
        mod     <- glm.fit(x, y, weights = w, family=gaussian(identity))
        # re-estimate variance function 
        mod.var <- glm.fit(z, mod$residuals^2, family=Gamma(log))

        dev0 <- dev1
        dev1 <- disp.dev(mod, mod.var)
    cat("Iteration ", i, ": deviance = ", format(dev1), "\n", sep="")

  }

  mod$call <- mean.formula
  mod.var$call <- var.formula
  class(mod) <- c("glm", "lm")
  class(mod.var) <- c("glm", "lm")
  result <- list(call=cl, mean=mod, var=mod.var, weights=w,
                 initial.deviance=initial.dev, deviance=dev1)
  class(result) <- "dispmod"
  return(result)

}

summary.dispmod <- function(object, ...)
{
  summary.mean <- summary(object$mean, dispersion=1, ...)
  summary.var  <- summary(object$var, dispersion=2, ...)
  cat("Normal dispersion model\n")
  cat(rep("-", options()$width), sep="")
  cat(paste(deparse(object$call), sep = "\n", 
                collapse = "\n"), "\n\n", sep = "")
  cat("Model for the mean \n")
  cat("------------------")
  print(summary.mean)
  cat("Model for the variance \n")
  cat("----------------------")
  print(summary.var)
  cat(paste("-2*logLik(max), constant var. =",
            format(object$initial.deviance),"\n"))
  cat(paste("-2*logLik(max), model         =", 
            format(object$deviance),"\n"))
  lrt <- object$initial.deviance - object$deviance
  df  <- object$var$df.null - object$var$df.residual
  cat(paste("LRT = ", format(lrt), " on ", df, " df, p-value = ", 
            format.pval(1-pchisq(lrt, df)), "\n", sep=""))        
  invisible(list(mean=summary.mean, var=summary.var, weights=object$weights,
                 initial.deviance=object$initial.deviance, 
                 deviance=object$deviance))          
}


##
##  Overdispersed binomial logit models
##

glm.binomial.disp <- function(object, maxit = 30, verbose = TRUE)
{
  if (class(object)[1] != "glm")
     stop("first argument must be a fitted model of class \"glm\" !")
  class <- class(object)
  if (!(family(object)$family == "binomial" & family(object)$link == "logit"))
     stop("overdispersed model fitting available only for \nbinomial regression models with logit link function!")

  pearson.X2 <- function(x) sum(residuals(x, "pearson")^2)
  
  y <- object$model[,1]       # observed proportion of success & failures
  trials <- apply(y, 1, sum)  # = object$prior.weights
  X <- model.matrix(object)
  p <- length(object$coefficients)
  n <- dim(X)[[1]]
  h <- lm.influence(object)$hat
  X2 <- pearson.X2(object)
  env <- parent.frame()
  assign("object", object, envir = env)
  
  # initial estimate of dispersion parameter
  phi <- max((X2 - (n-p)) / sum((trials-1)*(1-h)), 0)

  if(verbose)
    cat("\nBinomial overdispersed logit model fitting...\n")
  
  # loop until Pearson X2 approx equal to 1
  i <- 0
  converged <- TRUE
  while( abs(X2/(n-p)-1) > object$control$epsilon )
  {
    i <- i + 1
    if(i > maxit) 
      { converged <- FALSE
        break() }
    else if(verbose) 
      { cat("Iter. ", i, " phi:", format(phi), "\n") }

    # computes weights
    w <- 1/(1+phi*(trials-1))  
    # re-fit the model using update() evaluated in original model 
    assign("disp.weights", w, envir = env)
    object <- eval(expression(update(object, weights=disp.weights)), 
                   envir = env)
    #
    h <- lm.influence(object)$hat
    X2 <- pearson.X2(object)
    # current estimate of dispersion parameter
    phi <- max((X2 - sum(w*(1-h))) / sum(w*(trials-1)*(1-h)), 0)
  }

  if(verbose)
    { if(converged)
        { cat("Converged after", i, "iterations. \n")
          cat("Estimated dispersion parameter:", format(phi), "\n") 
          print(summary(object)) 
        }
      else
        warning("algoritm not converged after ", i, " iterations!")
    }

  object <- c(object, list(dispersion=phi, disp.weights=w))
  class(object) <- class
  invisible(object)
}

##
##  Overdispersed Poisson loglinear models
##

glm.poisson.disp <- function(object, maxit = 30, verbose = TRUE)
{
  if(class(object)[1] != "glm")
     stop("first argument must be a fitted model of class \"glm\" !")
  class <- class(object)
  if(!(family(object)$family == "poisson" & family(object)$link == "log"))
     stop("overdispersed model fitting available only for \npoisson regression models with log link function!")

  pearson.X2 <- function(x) sum(residuals(x, "pearson")^2)

  pw <- object$prior.weights
  y  <- object$y 
  mu <- object$fitted.values
  #R  <- object$R; Rinv <- solve(R); XtWXinv <- Rinv %*% t(Rinv)
  #X <- model.matrix(object); h <- diag( X %*% XtWXinv %*% t(X) )
  h <- lm.influence(object)$hat / object$weights
  p <- length(object$coefficients)
  n <- length(y)
  X2 <- pearson.X2(object)
  env <- parent.frame()
  assign("object", object, envir = env)

  # initial estimate of dispersion parameter
  phi <- max((X2 - (n-p)) / sum(mu*(1-mu*h)), 0)

  if(verbose)
    cat("\nPoisson overdispersed log-linear model fitting...\n")
    
  # loop until Pearson X2 approx equal to 1
  i <- 0
  converged <- TRUE
  while( abs(X2/(n-p)-1) > object$control$epsilon )
  {
    i <- i + 1
    if(i > maxit) 
      { converged <- FALSE
        break() }
    else if(verbose) 
      { cat("Iter. ", i, " phi:", format(phi), "\n") }

    # computes weights
    w <- 1/(1+(phi*mu))          
    # re-fit the model using update() evaluated in original model 
    assign("disp.weights", w, envir = env)
    object <- eval(expression(update(object, weights=disp.weights)), 
                   envir = env)

    mu <- object$fitted.values
    X2 <- pearson.X2(object)
    # current estimate of dispersion parameter
    phi <- max(sum( (y-mu)^2 / (mu*(mu+1/phi)) ) / (n-p), 0)
  }

  if(verbose)
    { if(converged)
        { cat("Converged after", i, "iterations. \n")
          cat("Estimated dispersion parameter:", format(phi), "\n") 
          print(summary(object)) 
        }
      else
        warning("algoritm not converged after ", i, " iterations!")
    }

  object <- c(object, list(dispersion=phi, disp.weights=w))
  class(object) <- class
  invisible(object)
}
