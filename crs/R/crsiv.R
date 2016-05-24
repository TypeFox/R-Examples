## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## alpha.min: minimum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## alpha.max: maximum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## ... optional arguments for crs()

## This function returns a list with some of the following elements:

## phi: the IV estimator of phi(y)
## alpha:  the Tikhonov regularization parameter
## phi.mat: the matrix with colums phi_1, phi_2 etc. over all iterations
## num.iterations:  the number of Landweber-Fridman iterations
## norm.stop: the vector of values of the objective function used for stopping
## norm.value: the norm not multiplied by the number of iterations
## convergence: a character string indicating whether/why iteration terminated

crsiv <- function(y,
                  z,
                  w,
                  x=NULL,
                  zeval=NULL,
                  weval=NULL,
                  xeval=NULL,
                  alpha=NULL,
                  alpha.min=1.0e-10,
                  alpha.max=1.0e-01,
                  alpha.tol=.Machine$double.eps^0.25,
                  deriv=0,
                  iterate.max=1000,
                  iterate.diff.tol=1.0e-08,
                  constant=0.5,
                  penalize.iteration=TRUE,
                  smooth.residuals=TRUE,
                  start.from=c("Eyz","EEywz"),
                  starting.values=NULL,
                  stop.on.increase=TRUE,
                  method=c("Landweber-Fridman","Tikhonov"),
                  opts=list("MAX_BB_EVAL"=10000,
                            "EPSILON"=.Machine$double.eps,
                            "INITIAL_MESH_SIZE"="r1.0e-01",
                            "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                            "MIN_POLL_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                            "DISPLAY_DEGREE"=0),
                  ...) {

  crs.messages <- getOption("crs.messages")

  ## This function was constructed initially by Samuele Centorrino
  ## <samuele.centorrino@univ-tlse1.fr>
  ## the following papers:

  ## A) Econometrica (2011) "Nonparametric Instrumental Regression"
  ## S. Darolles, Y. Fan, J.P. Florens, E. Renault, Volume 79,
  ## 1541-1565.

  ## B) Econometrics Journal (2010), volume 13, pp. S1–S27. doi:
  ## 10.1111/j.1368-423X.2010.00314.x "The practice of non-parametric
  ## estimation by solving inverse problems: the example of
  ## transformation models" Frederique Feve and Jean-Pierre Florens,
  ## IDEI and Toulouse School of Economics, Universite de Toulouse
  ## Capitole 21 alle de de Brienne, 31000 Toulouse, France. E-mails:
  ## feve@cict.fr, florens@cict.fr

  ## It was modified by Jeffrey S. Racine <racinej@mcmaster.ca> and all
  ## errors remain my responsibility. I am indebted to Samuele and the
  ## Toulouse School of Economics for their generous hospitality.

  ## First we require two functions, the first that conducts Regularized
  ## Tikhonov Regression' (aka Ridge Regression)

  ## This function conducts regularized Tikhonov regression which
  ## corresponds to (3.9) in Feve & Florens (2010).

  ## This function accepts as arguments

  ## alpha: penalty
  ## CZ:    row-normalized kernel weights for the `independent' variable
  ## CY:    row-normalized kernel weights for the `dependent' variable
  ## Cr:    row-normalized kernel weights for the `instrument/endogenous' variable (see NOTE below)
  ## r:     vector of conditional expectations (z can be E(Z|z) - see NOTE below)

  ## NOTE: for Cr, in the transformation model case treated in Feve &
  ## Florens (2010) this maps Z onto the Y space. In the IV case
  ## (Darrolles, Fan, Florens & Renault (2011) it maps W (the
  ## instrument) onto the space of the endogenous regressor Z.

  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.

  ## This function returns TBA (need better error checking!)

  ## phi:   the vector of estimated values for the unknown function at the evaluation points

  tikh <- function(alpha,CZ,CY,Cr.r){
    return(chol2inv(chol(alpha*diag(length(Cr.r)) + CY%*%CZ) %*% Cr.r)) ## This must be computable via ridge... step 1, step 2, same alpha...
  }

  ## This function applies the iterated Tikhonov approach which
  ## corresponds to (3.10) in Feve & Florens (2010).

  ## This function accepts as arguments

  ## alpha: penalty
  ## CZ:    row-normalized kernel weights for the `independent' variable
  ## CY:    row-normalized kernel weights for the `dependent' variable
  ## Cr:    row-normalized kernel weights for the `instrument/endogenous' variable (see NOTE below)
  ## r:     vector of conditional expectations (z can be E(Z|z) - see NOTE below)

  ## NOTE: for Cr, in the transformation model case treated in Feve &
  ## Florens (2010) this maps Z onto the Y space. In the IV case
  ## (Darrolles, Fan, Florens & Renault (2011) it maps W (the
  ## instrument) onto the space of the endogenous regressor Z.

  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.

  ## This function returns TBA (need better error checking!)

  ## phi:   the vector of estimated values for the unknown function at the evaluation points

  ## SSalpha: (scalar) value of the sum of square residuals criterion
  ## which is a function of alpha (see (3.10) of Feve & Florens (2010)

  ## Cr.r is always E.E.y.w.z, r is always E.y.w

  ittik <- function(alpha,CZ,CY,Cr.r,r) {
    invmat <- chol2inv(chol(alpha*diag(length(Cr.r)) + CY%*%CZ))
    tikh.val <- invmat %*% Cr.r
    phi <- tikh.val + alpha * invmat %*% tikh.val ## Not sure about this...
    return(sum((CZ%*%phi - r)^2)/alpha)     ## This is a sum of squared values so CZ%*%phi can be computed with fitted(crs())...
  }

  console <- newLineConsole()

  ## Basic error checking

  start.from <- match.arg(start.from)
  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")
  if(constant <= 0 || constant >=1) stop("constant must lie in (0,1)")
  if(iterate.diff.tol < 0) stop("iterate.diff.tol must be non-negative")

  ## Cast as data frames

  w <- data.frame(w)
  z <- data.frame(z)
  if(!is.null(x)) x <- data.frame(x)

  ## Check for evaluation data

  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w
  if(!is.null(x) && is.null(xeval)) xeval <- x

  method <- match.arg(method)

  if(!is.null(alpha) && alpha <= 0) stop("alpha must be positive")

  ## Set up formulas for multivariate w, z, and x if provided

  wnames <- names(w)
  znames <- names(z)
  names(weval) <- wnames
  names(zeval) <- znames

  ## If there exist exogenous regressors X, append these to the
  ## formulas involving Z (can be manually added to W by the user if
  ## desired)

  if(!is.null(x)) {
    xnames <- names(x)
    names(xeval) <- xnames
  }

  ## Now create evaluation data

  if(is.null(x)) {
    traindata <- data.frame(y,z,w)
    evaldata <- data.frame(zeval,weval)
  } else {
    traindata <- data.frame(y,z,w,x)
    evaldata <- data.frame(zeval,weval,xeval)
  }

  formula.yw <- as.formula(paste("y ~ ", paste(wnames, collapse= "+")))
  formula.phiw <- as.formula(paste("phi ~ ", paste(wnames, collapse= "+")))
  formula.residw <- as.formula(paste("(y-phi) ~ ", paste(wnames, collapse= "+")))

  if(is.null(x)) {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+")))
    formula.Ephiwz <- as.formula(paste("E.phi.w ~ ", paste(znames, collapse= "+")))
    formula.residwz <- as.formula(paste("residw ~ ", paste(znames, collapse= "+")))
  } else {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Ephiwz <- as.formula(paste("E.phi.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.residwz <- as.formula(paste("residw ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
  }

  if(!is.null(starting.values) && (NROW(starting.values) != NROW(evaldata))) stop(paste("starting.values must be of length",NROW(evaldata)))

  if(method=="Tikhonov") {

    ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
    ## bandwidths, one for y on w and one for phi(z) on w (in the
    ## first step we use E(y|w) as a proxy for phi(z) and use
    ## bandwidths for y on w).

    ## convergence flag returned for Landweber-Fridman, not Tikhonov,
    ## but value is required

    convergence <- NULL

    ## First we conduct the regression spline estimator of y on w

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weights and optimal smoothing for E(y|w)...", console)
    if(crs.messages) options(crs.messages=FALSE)
    model<-crs(formula.yw,opts=opts,data=traindata,...)
    if(crs.messages) options(crs.messages=TRUE)
    E.y.w <- predict(model,newdata=evaldata,...)
    B <- model.matrix(model$model.lm)
    KYW <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Next, we conduct the regression spline of E(y|w) on z

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing weights and optimal smoothing for E(E(y|w)|z)...", console)
    } else {
      console <- printPush("Computing weights and optimal smoothing for E(E(y|w)|z,x)...", console)
    }
    if(crs.messages) options(crs.messages=FALSE)
    model <- crs(formula.Eywz,opts=opts,data=traindata,...)
    if(crs.messages) options(crs.messages=TRUE)
    E.E.y.w.z <- predict(model,newdata=evaldata,...)
    B <- model.matrix(model$model.lm)
    KYWZ <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Next, we minimize the function ittik to obtain the optimal value
    ## of alpha (here we use the iterated Tikhonov function) to
    ## determine the optimal alpha for the non-iterated scheme. Note
    ## that the function `optimize' accepts bounds on the search (in
    ## this case alpha.min to alpha.max))

    ## E(r|z)=E(E(phi(z)|w)|z)
    ## \phi^\alpha = (\alpha I+CzCw)^{-1}Cr x r

    if(is.null(alpha)) {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush("Numerically solving for alpha...", console)
      alpha <- optimize(ittik, c(alpha.min,alpha.max), tol = alpha.tol, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    }

    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha to get a first stage estimate of phi

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing initial phi(z) estimate...", console)
    } else {
      console <- printPush("Computing initial phi(z,x) estimate...", console)
    }
    phi <- as.vector(tikh(alpha, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z))

    ## KYWZ and KZWS no longer used, save memory

    rm(KYW,KYWZ)

    ## Conduct kernel regression of phi(z) on w

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing optimal smoothing and weights for E(phi(z)|w)...", console)
    } else {
      console <- printPush("Computing optimal smoothing and weights for E(phi(z,x)|w)...", console)
    }
    if(crs.messages) options(crs.messages=FALSE)
    model <- crs(formula.phiw,opts=opts,data=traindata,...)
    if(crs.messages) options(crs.messages=TRUE)
    E.phi.w <- predict(model,newdata=evaldata,...)
    B <- model.matrix(model$model.lm)
    KPHIW <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Conduct kernel regression of E(phi(z)|w) on z

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing optimal smoothing and weights for E(E(phi(z)|w)|z)...", console)
    } else {
      console <- printPush("Computing optimal smoothing and weights for E(E(phi(z,x)|w)|z,x)...", console)
    }
    if(crs.messages) options(crs.messages=FALSE)
    model <- crs(formula.Ephiwz,opts=opts,data=traindata,...)
    if(crs.messages) options(crs.messages=TRUE)
    B <- model.matrix(model$model.lm)
    KPHIWZ <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Next, we minimize the function ittik to obtain the optimal value of
    ## alpha (here we use the iterated Tikhonov approach) to determine the
    ## optimal alpha for the non-iterated scheme.

    if(is.null(alpha)) {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush("Iterating and computing the numerical solution for alpha...", console)
      alpha <- optimize(ittik,c(alpha.min,alpha.max), tol = alpha.tol, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    }

    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha.

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing final phi(z) estimate...", console)
    } else {
      console <- printPush("Computing final phi(z,x) estimate...", console)
    }
    phi <- as.vector(tikh(alpha, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z))

    console <- printClear(console)
    console <- printPop(console)

    if((alpha-alpha.min)/alpha.min < 0.01) warning(paste(" Tikhonov parameter alpha (",formatC(alpha,digits=4,format="f"),") is close to the search minimum (",alpha.min,")",sep=""))
    if((alpha.max-alpha)/alpha.max < 0.01) warning(paste(" Tikhonov parameter alpha (",formatC(alpha,digits=4,format="f"),") is close to the search maximum (",alpha.max,")",sep=""))

    ## phi.0 is the conditional mean model. We compute lambda =
    ## fitted(phi.0)-phi then transform y via
    ## y.lambda=y-lambda. Here we overwrite y so that we can reuse the
    ## formula. Before that, save the proper residuals and then push
    ## these into the model. June 9 2011 - I am concerned because
    ## phi and the fitted values from this approach are _identical_
    ## (I expected approximately equal).

    ## Feb 21 2012 - JP Florens said the starting point should be
    ## E[E[Y|W]|Z], below we do E[Y|Z]... certainly works, could we
    ## shorten the iterative process?

    if(crs.messages) options(crs.messages=FALSE)
    phi.0 <- crs(formula.yz,opts=opts,data=traindata,...)

    residuals.phi <- traindata$y-phi
    traindata$y <- traindata$y - (fitted(phi.0)-phi)

    model <- crs(formula.yz,
                 cv="none",
                 degree=phi.0$degree,
                 segments=phi.0$segments,
                 lambda=phi.0$lambda,
                 include=phi.0$include,
                 kernel=phi.0$kernel,
                 basis=phi.0$basis,
                 knots=phi.0$knots,
                 tau=phi.0$tau,
                 deriv=deriv,
                 data=traindata)
    if(crs.messages) options(crs.messages=TRUE)

    model$residuals <- residuals.phi
    model$phi <- phi
    model$alpha <- alpha

    return(model)

  } else {

    ## Landweber-Fridman

    ## Create storage vector/matrix

    norm.stop <- numeric()

    ## Compute E(Y|w) for the stopping rule

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing and E(Y|w) for the stopping rule...",sep=""),console)

    if(crs.messages) options(crs.messages=FALSE)
    model.E.y.w <- crs(formula.yw,opts=opts,data=traindata,...)
    E.y.w <- predict(model.E.y.w,newdata=evaldata,...)
    if(crs.messages) options(crs.messages=TRUE)

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing and phi(z,x) for iteration 1...",sep=""),console)
    }

    ## Initial value taken from E(E(Y|w)|z) or E(Y|z) or overridden
    ## and passed in, formulae all operate on phi. phi.0.NULL flag set

    if(crs.messages) options(crs.messages=FALSE)
    if(is.null(starting.values)) {
      phi.0.NULL <- TRUE
      phi.0 <- crs(formula.yz,opts=opts,data=traindata,...)
      ## First compute phi.0 (not passed in) then phi
      if(start.from == "Eyz") {
        ## Start from E(Y|z)
        phi <- predict(phi.0,newdata=evaldata,...)
      } else {
        ## Start from E(E(Y|w)|z)
        E.y.w <- fitted(crs(formula.yw,opts=opts,data=traindata,...))
        model.E.E.y.w.z <- crs(formula.Eywz,opts=opts,data=traindata,...)
        phi <- predict(model.E.E.y.w.z,newdata=evaldata,...)
      }
    } else {
      phi.0.NULL <- FALSE
      phi.0.input <- starting.values
      ## First compute phi (passed in) then phi.0
      phi <- starting.values
      phi.0 <- crs(formula.yz,opts=opts,data=traindata,...)
    }

    starting.values.phi <- phi
    
    if(crs.messages) options(crs.messages=TRUE)

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing for E(Y-phi(z)|w) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing  for E(Y-phi(z,x)|w) for iteration 1...",sep=""),console)
    }
    if(crs.messages) options(crs.messages=FALSE)
    if(smooth.residuals) {
      model.residw <- crs(formula.residw,opts=opts,data=traindata,...)
      residw <- predict(model.residw,newdata=evaldata,...)
      model.predict.residw.z <- crs(formula.residwz,opts=opts,data=traindata,...)
    } else {
      model.E.phi.w <- crs(formula.phiw,opts=opts,data=traindata,...)
      residw <- predict(model.E.y.w,newdata=evaldata,...)-predict(model.E.phi.w,newdata=evaldata,...)
      model.predict.residw.z <- crs(formula.residwz,opts=opts,data=traindata,...)
    }
    if(crs.messages) options(crs.messages=TRUE)

    if(phi.0.NULL) {
      phi <- predict(phi.0,newdata=evaldata,...) + constant*predict(model.predict.residw.z,newdata=evaldata,...)
    } else {
      phi <- phi.0.input + constant*predict(model.predict.residw.z,newdata=evaldata,...)
    }

    phi.mat <- phi
    norm.stop[1] <- sum(residw^2)/sum(E.y.w^2)

    for(j in 2:iterate.max) {

      console <- printClear(console)
      console <- printPop(console)
      if(is.null(x)) {
        console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration ", j,"...",sep=""),console)
      } else {
        console <- printPush(paste("Computing optimal smoothing and phi(z,x) for iteration ", j,"...",sep=""),console)
      }

      if(crs.messages) options(crs.messages=FALSE)
      if(smooth.residuals) {
        model.residw <- crs(formula.residw,opts=opts,data=traindata,...)
        residw <- predict(model.residw,newdata=evaldata,...)
        model.predict.residw.z <- crs(formula.residwz,opts=opts,data=traindata,...)
      } else {
        model.E.phi.w <- crs(formula.phiw,opts=opts,data=traindata,...)
        residw <- predict(model.E.y.w,newdata=evaldata,...)-predict(model.E.phi.w,newdata=evaldata,...)
        model.predict.residw.z <- crs(formula.residwz,opts=opts,data=traindata,...)
      }
      if(crs.messages) options(crs.messages=TRUE)

      phi <- phi + constant*predict(model.predict.residw.z,newdata=evaldata,...)
      phi.mat <- cbind(phi.mat,phi)

      norm.stop[j] <- ifelse(penalize.iteration,j*sum(residw^2)/sum(E.y.w^2),sum(residw^2)/sum(E.y.w^2))

      ## The number of iterations in LF is asymptotically equivalent
      ## to 1/alpha (where alpha is the regularization parameter in
      ## Tikhonov).  Plus the criterion function we use is increasing
      ## for very small number of iterations. So we need a threshold
      ## after which we can pretty much confidently say that the
      ## stopping criterion is decreasing.  In Darolles et al. (2011)
      ## \alpha ~ O(N^(-1/(min(beta,2)+2)), where beta is the so
      ## called qualification of your regularization method. Take the
      ## worst case in which beta = 0 and then the number of
      ## iterations is ~ N^0.5.

      if(j > round(sqrt(nrow(traindata))) && !is.monotone.increasing(norm.stop)) {

        ## If stopping rule criterion increases or we are below stopping
        ## tolerance then break

        if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
          convergence <- "STOP_ON_INCREASE"
          break()
        }
        if(abs(norm.stop[j-1]-norm.stop[j]) < iterate.diff.tol) {
          convergence <- "ITERATE_DIFF_TOL"
          break()
        }

      }

      convergence <- "ITERATE_MAX"

    }

    norm.value <- norm.stop/(1:length(norm.stop))

    ## Extract minimum, and check for monotone increasing function and
    ## issue warning in that case. Otherwise allow for an increasing
    ## then decreasing (and potentially increasing thereafter) portion
    ## of the stopping function, ignore the initial increasing portion,
    ## and take the min from where the initial inflection point occurs
    ## to the length of norm.stop

    if(which.min(norm.stop) == 1 && is.monotone.increasing(norm.stop)) {
      warning("Stopping rule increases monotonically (consult model$norm.stop):\nThis could be the result of an inspired initial value (unlikely)\nNote: we suggest manually choosing phi.0 and restarting (e.g. instead set `starting.values' to E[E(Y|w)|z])")
      convergence <- "FAILURE_MONOTONE_INCREASING"
      ## Ignore the initial increasing portion, take the min to the
      ## right of where the initial inflection point occurs
      j <- 1
      while(norm.value[j+1] > norm.value[j]) j <- j + 1
      j <- j-1 + which.min(norm.value[j:length(norm.value)])
      phi <- phi.mat[,j]
#      phi <- starting.values.phi
    } else {
      ## Ignore the initial increasing portion, take the min to the
      ## right of where the initial inflection point occurs
      j <- 1
      while(norm.stop[j+1] > norm.stop[j]) j <- j + 1
      j <- j-1 + which.min(norm.stop[j:length(norm.stop)])
      phi <- phi.mat[,j]
    }

    ## phi.0 is the conditional mean model. We compute lambda =
    ## fitted(phi.0)-phi then transform y via
    ## y.lambda=y-lambda. Here we overwrite y so that we can reuse the
    ## formula. Before that, save the proper residuals and then push
    ## these into the model. June 9 2011 - I am concerned because
    ## phi and the fitted values from this approach are _identical_
    ## (I expected approximately equal).

    residuals.phi <- traindata$y-phi
    traindata$y <- traindata$y - (fitted(phi.0)-phi)

    if(crs.messages) options(crs.messages=FALSE)
    model <- crs(formula.yz,
                 cv="none",
                 degree=phi.0$degree,
                 segments=phi.0$segments,
                 lambda=phi.0$lambda,
                 include=phi.0$include,
                 kernel=phi.0$kernel,
                 basis=phi.0$basis,
                 knots=phi.0$knots,
                 tau=phi.0$tau,
                 deriv=deriv,
                 data=traindata)
    if(crs.messages) options(crs.messages=TRUE)

    model$residuals <- residuals.phi
    model$phi <- phi
    model$phi.mat <- phi.mat
    model$num.iterations <- j
    model$norm.stop <- norm.stop
    model$norm.value <- norm.value
    model$convergence <- convergence
    model$starting.values.phi <- starting.values.phi
    console <- printClear(console)
    console <- printPop(console)

    if(j == iterate.max) warning(" iterate.max reached: increase iterate.max or inspect norm.stop vector")

    return(model)

  }

}
