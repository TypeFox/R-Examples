########## R-function: asp2 ##########

# For fitting adaptive semiparametric regression models.

"asp2" <-   function (form, #adap = FALSE, #random = NULL, group = NULL, family = "gaussian",
            spar.method = "REML", contrasts=NULL, omit.missing = NULL,returnFit=FALSE, niter = 20, niter.var = 50,tol=1e-6, tol.theta=1e-6,
          #  weights = NULL,correlation=NULL,
            control=NULL){

  epsilon.fit=epsilon.theta=NULL


#if (family != "gaussian") stop("Only Gaussian response supported. Use AdaptFit for GAMs (without simultaneous confidence bands).")
#if (!is.null(correlation)) stop("Correlated errors not supported. Use AdaptFit for correlated errors (without simultaneous confidence bands).")
#if (!is.null(random)) stop("Random effects not supported. Use AdaptFit (without simultaneous confidence bands).")

#read the formula

  spm.info <- aspmFormReadOS(form, omit.missing,constrasts=contrasts)
  design.info <- spmDesignOS(spm.info)

#get the design matrices

  Xb <- design.info$X
  Zb <- design.info$Z
  Wb <- cbind(Xb, Zb)
  asp.info <- NULL
  model.matrices <- list(Xb = Xb, Zb = Zb, Wb = Wb)

#define weights for grouped binomial data
#  if ((family == "binomial") & (is.null(weights)))
#    weights = rep(1, nrow(Xb))
#  if (!is.null(dim(spm.info$y))) {
#    weights <- apply(spm.info$y, 1, sum)
#    spm.info$y <- (spm.info$y)[, 1]/weights
#  }

#get non-adaptive fit 
  aspm.obj <- aspmOS(spm.info, random=NULL, group=NULL, family="gaussian", spar.method, omit.missing,correlation=NULL,weights=NULL,control=control)

#  if (family != "gaussian") aspm.obj$fit$fitted <- NULL
  fit <- aspm.obj$fit$fitted

#redefine correlation matrix in case of correlated errors

#  if(!is.null(correlation))     y.cov <- aspm.obj$fit$sigma^2*corMatrix(aspm.obj$fit$modelStruct$corStruct)
#  else
  if(!is.null(aspm.obj$fit$sigma))   y.cov <- aspm.obj$fit$sigma^2#*diag(rep(1,nrow(Xb)))
  else y.cov <- NULL

#get the initial estimates for the variance of random effects
  b.hat <- as.vector(unlist(aspm.obj$fit$coef$random))
  kc <- kb <- theta <- sigma.theta <- knots.theta <- adap.ind <- theta.ind <- NULL
  Wc <- list()
  stb <- 1
  l <- 0
  adap=FALSE
  if (!is.null(spm.info$pen))     adap.ind <- spm.info$pen$adap
  if (!is.null(spm.info$krige))   adap.ind <- c(adap.ind, spm.info$krige$adap)
  if (sum(adap.ind) >= 1)       adap=TRUE
  if ((!is.null(spm.info$pen) | !is.null(spm.info$krige)) & (adap == TRUE)) {
    if (!is.null(spm.info$pen)) 
      for (l in 1:length(spm.info$pen$name)) {
        knot <- design.info$spm.info$pen$knots[[l]]
        var.knots <- spm.info$pen$var.knots[[l]]
        if (spm.info$pen$var.basis[[l]]== "trunc.poly"){cat("\nSetting var.basis=\"tps\",var.degree=3\n\n");spm.info$pen$var.basis[[l]]="tps";spm.info$pen$var.degree[[l]]=3}
        if (!is.null(var.knots))
            design.info.theta <- spmDesignOS(aspmFormReadOS(as.formula(paste("1~f(knot,basis='",spm.info$pen$var.basis[[l]],"',degree=",spm.info$pen$var.degree[[l]],",knots=var.knots)",sep="")), omit.missing))
        else design.info.theta <- spmDesignOS(aspmFormReadOS(as.formula(paste("1~f(knot,basis='",spm.info$pen$var.basis[[l]],"',degree=",spm.info$pen$var.degree[[l]],")",sep="")), omit.missing))
        Xc <- design.info.theta$X
        Zc <- design.info.theta$Z
        Wc[[l]] <- cbind(Xc, Zc)
        kc <- c(kc, ncol(design.info.theta$X), ncol(design.info.theta$Z))
        kb <- c(kb, nrow(design.info.theta$Z))
        all <- rep(1, nrow(Zc))
        bb <- b.hat[stb:(kb[l] + stb - 1)]
        b.hat.var <- log((bb - mean(bb))^2/var(bb))
        b.fit <- lme(b.hat.var ~ Xc - 1, random = list(all = pdIdent(~Zc - 1)))
        theta <- c(theta, c(as.vector(b.fit$coef$fixed), as.vector(unlist(b.fit$coef$random))))
        theta.ind <- c(theta.ind, rep(adap.ind[l], ncol(Wc[[l]])))
        sigma.theta <- c(sigma.theta,  as.vector(b.fit$sigma^2 * exp(2 * unlist(b.fit$modelStruct$reStruct))))
        knots.theta <- c(knots.theta, list(design.info.theta$spm.info$pen$knots[[1]]))
        stb <- stb + kb[l]
      }
    if (!is.null(spm.info$krige)) {
stop("Surface estimation currently not supported in AdaptFitOS. Use AdaptFit")
#      knot <- spm.info$krige$knots
#      var.knots <- spm.info$krige$var.knot
#      knots1 <- knot[, 1]
#      knots2 <- knot[, 2]
#      theta.frame <- list(knots1 = knots1, knots2 = knots2,
#                          var.knots = var.knots)
#      attach(theta.frame)
#      if (!is.null(var.knots))
#        design.info.theta <- spmDesignOS(aspmFormReadOS(as.formula(paste("1~f(knots1,knots2,knots=var.knots)")), omit.missing))
#      else design.info.theta <- spmDesignOS(aspmFormReadOS(as.formula(paste("1~f(knots1,knots2)")), omit.missing))
#      Xc <- design.info.theta$X
#      Zc <- design.info.theta$Z
#      Wc[[l + 1]] <- cbind(Xc, Zc)
#      kc <- c(kc, ncol(design.info.theta$X), ncol(design.info.theta$Z))
#      kb <- c(kb, nrow(design.info.theta$Z))
#      all <- rep(1, nrow(Zc))
#      assign("Xc", Xc)
#      assign("Zc", Zc)
#      bb <- b.hat[stb:(kb[l + 1] + stb - 1)]
#      b.hat.var <- log((bb - mean(bb))^2/var(bb))
#      b.fit <- lme(b.hat.var ~ Xc - 1, random = list(all = pdIdent(~Zc - 1)))
#      theta <- c(theta, c(as.vector(b.fit$coef$fixed), as.vector(unlist(b.fit$coef$random))))
#      theta.ind <- c(theta.ind, rep(adap.ind[l + 1], length(theta)))
#      sigma.theta <- c(sigma.theta, as.vector(b.fit$sigma^2 * exp(2 * unlist(b.fit$modelStruct$reStruct))))
#      knots.theta <- c(knots.theta, list(design.info.theta$spm.info$krige$knots[[1]]))
#      detach(theta.frame)
    }


#store the model matrices
    
    Wc <- bdiag(Wc)
#  if (family != "gaussian")
#    {
#      if (family == "poisson")
#      V <- c(fit)
#    if (family == "binomial")
#      V <- c(fit * (1 - fit) * weights)
#    V12 <- sqrt(V)
#    Zb <- t(t(Zb) * V12)
#    }
    ZVZ <- t(Zb)%*%Zb
    XVX <- t(Xb)%*%Xb
    ZVX <- t(Zb)%*%Xb
    
    model.matrices <- list(Xb = Xb, Zb = Zb, Wc = Wc, ZVZ=ZVZ, XVX=XVX, ZVX=ZVX, kb = kb, kc = kc)

#iterate for the next estimate of the variance of random effects

    for (i in 1:niter) {
      theta.info <- list(theta = theta, sigma.theta = sigma.theta, asp.info = aspm.obj, model.matrices = model.matrices)

	#get the new estimate for the spline coefficients

      theta.obj <- aspGetThetaOS(theta.info, niter.var, family="gaussian", weights=NULL,spar.method,returnFit,tol.theta=tol.theta)

	#obtain variance of random effects

      sigma.theta <- theta.obj$sigma.theta
      theta <- theta.obj$theta
      theta1 <- c(theta)
      theta1[!theta.ind] <- 0
      ran.var <- exp(c(Wc %*% theta1))
      
	#obtain remaining estimate with the estimated variance of random effects

      aspm.obj <- aspmOS(spm.info, random=NULL, group=NULL, family="gaussian", spar.method, omit.missing, Si.b = ran.var, weights = NULL,correlation=NULL,control=control)
#      if (family != "gaussian")  aspm.obj$fit$fitted <- NULL
      fit1 <- aspm.obj$fit$fitted
      epsilon.fit <- sum((fit - fit1)^2)/sum(fit^2)
      epsilon.theta=theta.obj$epsilon.theta
      fit <- fit1

  #exit if epsilon is small
      if (epsilon.fit <= tol) break
      if (i == niter)
          if (returnFit==TRUE) {cat("Warning: Iteration limit reached without convergence");break}
          else  stop ("Iteration limit reached without convergence")
    }

#extract objects for the output related to the adaptive fit

    theta <- theta1
    subknots.pen <- knots.theta
    subknots.krige <- design.info.theta$info$krige$knots
    random.var.fix <- as.vector(aspm.obj$fit$sigma^2 * exp(2 * unlist(aspm.obj$fit$modelStruct$reStruct)))
    coeff.random.var <- c(theta) + c(ginv(Wc) %*% rep(log(random.var.fix), kb))
    model.matrices <- list(Xb = Xb, Zb = Zb, Wb = cbind(Xb, Zb), Wc = Wc)
    #if(!is.null(correlation))       y.cov <- aspm.obj$fit$sigma^2*corMatrix(aspm.obj$fit$modelStruct$corStruct)
    #else
     y.cov <- aspm.obj$fit$sigma^2#*diag(rep(1,nrow(Xb)))
    asp.info <- list(subknots = list(subknots.pen = subknots.pen, subknots.krige = subknots.krige), coef.random.var = coeff.random.var,  var.random.var = sigma.theta)
  }

#extract remaining objects for the output

  fitted <- aspm.obj$fit$fitted
  coeff.mean.fixed <- aspm.obj$fit$coef$fixed
  coeff.mean.random <- as.vector(unlist(aspm.obj$fit$coef$random))
  lin.x <- aspm.obj$info$lin$x
  pen.x <- aspm.obj$info$pen$x
  krige.x <- aspm.obj$info$krige$x
  pen.knots <- aspm.obj$info$pen$knots
  krige.knots <- aspm.obj$info$krige$knots
  random.var <- diag(aspm.obj$aux$random.var)


if (adap==FALSE) aspm.obj$info$pen$adap=rep(FALSE,length(spm.info$pen$adap))


#construct ouput

  asp.info <- c(list(fitted = fitted, coef.mean = c(coeff.mean.fixed, coeff.mean.random), design.matrices = model.matrices, 
                     x = list(lin.x = lin.x, pen.x = pen.x, krige.x = krige.x), 
                     knots = list(pen.knots = pen.knots, krige.knots = krige.knots), 
                     y.cov=y.cov, random.var = random.var,epsilon.fit=epsilon.fit,epsilon.theta=epsilon.theta),asp.info)

  asp.fit.obj <- c(asp.info, aspm.obj)
  class(asp.fit.obj) <-"asp"
  return(asp.fit.obj)
}
