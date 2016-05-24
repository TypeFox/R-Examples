
KFconvar <- function(model,
  P0cov = FALSE, barrier = list(type = "1", mu = 0), debug = TRUE)
{
  #NOTE using 'rescale=FALSE' uses relative variances and 'cpar=1'
  #NOTE see use factor (n/(n-1)) * s2

  # compute the value of the concentrated parameter

  n <- length(model@y)
  ss <- char2numeric(model, P0cov = P0cov, rescale = FALSE)

  mod <- list(Z = ss$Z, a = ss$a0, P = ss$P0, T = ss$T, 
    V = ss$Q, h = ss$H, Pn = ss$P0)
  res <- stats::KalmanLike(model@y, mod, -1, TRUE)
  s2 <- as.vector(res[[2]])
  mll <- 0.5 * n * log(2 * pi + 1) + n * as.vector(res[[1]])

  if (barrier$mu != 0)
  {
    bar <- barrier.eval(model, barrier$type, barrier$mu, 
      FALSE, FALSE)$barrier
  } else bar <- 0

  mll <- mll + bar

  if (debug && barrier$mu == 0)
  {
##NOTE
#do not update KF for new value of s2 or new variances, 
#otherwise s2 depends on the remaining variances

    kf <- KF(model@y, ss) #convergence, t0
    s2b <- sum(kf$v^2 / kf$f) / n #length(model@diffy)
    nh <- n / 2
    mllb <- nh * log(2 * pi + 1) + 0.5 * sum(log(kf$f)) + nh * log(s2b)

    stopifnot(all.equal(s2, s2b))
    stopifnot(all.equal(mll, mllb))
  }

  list(mll = mll, cpar = s2)
}

mloglik.td <- function(x, model, 
  #KF.version = c("KFKSDS", "stats", "StructTS", "KFAS", "FKF", "sspir", "dlm", "dse"),
  KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version),
  KF.args = list(), check.KF.args = TRUE,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999)
{
  if (!missing(x)) if(!is.null(x))
    model@pars <- x

  if (!is.null(model@xreg)) {
    y2 <- model@y - model@xreg %*% cbind(model@pars[model@ss$xreg])
  } else 
    y2 <- model@y

  P0cov <- if (is.null(KF.args$P0cov)) FALSE else KF.args$P0cov

  if (!is.null(model@cpar))
  {
    mll <- KFconvar(model, P0cov = P0cov, barrier = barrier, debug = TRUE)$mll

    # "KF.version" and "check.KF.args" are ignored when a parameter 
    # is concentrated out of likelihood function
    if (KF.version != "KFKSDS")
      warning("argument 'KF.version = ", KF.version, "' was ignored.")

    return(mll)
  }

  # Kalman filter
  # evaluate the minus log-likelihood function

  mll <- KalmanFilter(y = y2, ss = char2numeric(model, P0cov),
    KF.version = KF.version, KF.args = KF.args,
    check.args = check.KF.args, debug = FALSE)$mloglik

  # barrier term (optional)

  if (barrier$mu != 0)
  {
    bar <- barrier.eval(model, barrier$type, barrier$mu, 
      FALSE, FALSE)$barrier
  } else bar <- 0

  mll <- mll + bar

  if (!is.finite(mll))
    mll <- sign(mll) * inf

  # this is the minus log-lik, penalize with a large positive value
  if (is.na(mll))
    mll <- abs(inf)

  mll
}

mloglik.td.deriv <- function(model, gradient = TRUE, infomat = TRUE,
  KF.args = list(), version = c("1", "2"), kfres = NULL,
  convergence = c(0.001, length(model@y)))
{
  version <- match.arg(version)[1]

  if (!is.null(model@xreg))
  {
##FIXME TODO
    if (!is.null(model@transPars) && model@transPars != "square")
      stop("Analytical derivatives are not available for model@transPars=", 
        dQuote(model@transPars))

    xregnms <- model@ss$xreg
    idxreg <- match(xregnms, names(model@pars))
    varnms <- names(model@pars[-idxreg])
    #xreg <- list(xreg = xreg, coefs = model@pars[idxreg])
    # adjust model@y for "KF.deriv.C" and "transPars" 
    # "transPars" uses model@y only if model@transPars is "StructTS"
    model@y <- model@y - model@xreg %*% cbind(model@pars[idxreg])
  } else {
    varnms <- names(model@pars)
  }

  if (!is.null(model@transPars) && any(gradient, infomat))
  {
    #if (!is.null(xreg) && model@transPars == "StructTS")
    #  model@y <- model@y - xreg$xreg %*% cbind(xreg$coefs)
    dtrans <- transPars(model, gradient = TRUE, hessian = infomat)
    if (!is.null(model@xreg))
      dtrans$gradient[xregnms] <- 0
  } #else dtrans <- NULL

  if (is.null(kfres))
  {
    P0cov <- if (is.null(KF.args$P0cov)) FALSE else KF.args$P0cov
    ss <- char2numeric(model, P0cov)
    kf <- KF.deriv.C(model@y, ss, xreg = model@xreg, convergence = convergence)    
  } else { kf <- kfres }

  vof <- kf$v / kf$f
  invf <- 1 / kf$f

  # gradient

  g <- model@pars
  g[] <- NA

  if (gradient)
  {
    if (version == "1")
    {
      #for(i in seq(along = g))
      for(i in varnms)
      {
        part1 <- invf * kf$df[,i] * (1 - vof * kf$v)
        part2 <- kf$dv[,i] * vof
        g[i] <- sum(0.5 * part1 + part2)
      }
    } else if (version == "2")
    {
      #for(i in seq(along = g))
      for(i in varnms)
      {
        part1 <- kf$df[,i] / kf$f
        part2 <- kf$dv[,i] * vof
        part3 <- vof^2 * kf$df[,i]
        part4 <- vof * kf$dv[,i]
        g[i] <- 0.5 * sum(part1 + part2 - part3 + part4)
      }
    } else
      stop(paste("version", sQuote(version), 
        "is not implemented in 'mloglik.td.deriv'."))

    if (!is.null(model@transPars))
    {
      #g[varnms] <- g[varnms] * dtrans$gradient
      # dtrans$gradient[xregnms] is set to zeros above
      g <- g * dtrans$gradient
    }

    if (!is.null(model@xreg))
    {
      if (length(xregnms) == 1) {
        g[xregnms] <- sum(kf$dv[,xregnms] * vof)
      } else
        g[xregnms] <- colSums(kf$dv[,xregnms] * vof)
    }
  }

  # information matrix

  IM <- matrix(nrow = length(g), ncol = length(g))
  rownames(IM) <- colnames(IM) <- names(g)

  if (infomat)
  {
    invfsq <- invf^2
    #for(i in seq(np)) for(j in seq(np))
    #for (i in varnms) for (j in varnms)
    tmp <- cbind(rbind(varnms, varnms), combn(varnms,2))
    for (ij in seq.int(ncol(tmp)))
    {
      i <- tmp[1,ij]
      j <- tmp[2,ij]
      IM[i,j] <- sum(0.5 * invfsq * kf$df[,i] * kf$df[,j])
      IM[i,j] <- IM[i,j] + sum(kf$dv[,i] * kf$dv[,j] * invf)
      if (i != j)
        IM[j,i] <- IM[i,j]
    }

    if (!is.null(model@transPars))
    {
      #IM[varnms,varnms] <- tcrossprod(dtrans$gradient) * IM[varnms,varnms]
      IM <- tcrossprod(dtrans$gradient) * IM
    }

    if (!is.null(model@xreg))
    {
      allnms <- rownames(IM)
      for (i in allnms) for (j in xregnms)
        IM[j,i] <- IM[i,j] <- sum(kf$dv[,i] * kf$dv[,j] * invf)
    }

    #isSymmetric(IM)
  }

  list(gradient = g, infomat = IM)
}

mloglik.td.grad <- function(x, model, KF.version, KF.args = list(), 
  convergence = c(0.001, length(model@y)),
  check.KF.args, barrier, inf)
{
  if (!missing(x)) if(!is.null(x))
    model@pars[] <- x

  g <- mloglik.td.deriv(model = model, gradient = TRUE, infomat = FALSE,
    KF.args = KF.args, version = "1", kfres = NULL,
    convergence = convergence)$gradient

  isna <- is.na(g)
  if (any(isna))
    g[isna] <- inf

  g
}
