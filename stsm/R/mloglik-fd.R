
##NOTE
# argument "xreg" is not used in "mloglik.fd"
# this argument is required when used within "optim"
# "maxlik.fd.optim" where "xreg" is passed to "mloglik.fd.deriv"
# when analytical derivatives are used; in that case "xreg" is a list
# of constant terms, it is not the matrix of external regressors

mloglik.fd <- function(x, model, 
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, xreg)
{
  if (!missing(x)) if (!is.null(x))
    model@pars <- x

  n <- length(model@diffy)
  pi2 <- 2 * pi

  if (!is.null(model@xreg))
  {
    y2 <- model@y - model@xreg %*% cbind(model@pars[model@ss$xreg])
    model@diffy <- model@fdiff(y2, frequency(y2))
    # in general update only if (!is.null(model@ssd))
    model@ssd <- as.vector(Mod(fft(model@diffy))^2 / (pi2 * n))
  }

  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (pi2 * n)
  } else
    pg <- model@ssd

  # spectral generating function of the stationary differenced data

  sgf <- stsm.sgf(model, FALSE, FALSE, FALSE)$sgf

  # minus log-likelihood

  if (sgf[1] == 0)
  {
    sgf <- sgf[-1]
    pg <- pg[-1]
  }

  mll <- 0.5 * (n * log(pi2) + sum(log(sgf))) + pi * sum(pg/sgf)

  # barrier term (optional)

  if (barrier$mu != 0)
  {
    bar <- barrier.eval(model, barrier$type, barrier$mu, 
      FALSE, FALSE)$barrier
  } else bar <- 0

  mll <- mll + bar

  if (!is.finite(mll))
    mll <- sign(mll) * inf

  if (is.na(mll))
    mll <- abs(inf) #minus log-lik, penalize with a large value

  mll
}

mloglik.fd.deriv <- function(model, xreg = NULL,
  gradient = TRUE, hessian = TRUE, infomat = TRUE, modcovgrad = TRUE,
  barrier = list(type = c("1", "2"), mu = 0),
  version = c("2", "1"))
{
  version <- match.arg(version)[1]

  dtpars <- if (version == "1") TRUE else FALSE
  if (!is.null(model@transPars) && version == "2" && any(gradient, hessian, modcovgrad)) {
    dtrans <- transPars(model, gradient = gradient, hessian = hessian)
  } else dtrans <- NULL

  pi2 <- 2 * pi

  tmp <- stsm.sgf(model, TRUE, (hessian || modcovgrad), dtpars)
  sgf <- tmp$sgf
  sgf.d1 <- tmp$gradient
  sgf.d2 <- tmp$hessian

  if (!is.null(model@xreg))
  {
    n <- length(model@diffy)
    xregnms <- model@ss$xreg
    xregcoefs <- model@pars[xregnms]
    ncxreg <- ncol(model@xreg)

    #y2 <- model@y - xregmat %*% cbind(xregcoefs)
    y2 <- model@y - model@xreg %*% cbind(xregcoefs)
    tmp <- fft(model@fdiff(y2, frequency(y2)))
    model@ssd <- as.numeric(Mod(tmp)^2 / (pi2 * n))

  } else ncxreg <- 0

  if (any(c(gradient, hessian, infomat, modcovgrad)) && barrier$mu != 0)
    bar <- barrier.eval(model, barrier$type, barrier$mu, 
      gradient, (hessian || infomat))

  # if "xreg" is not NULL model@ssd is already defined above
  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (pi2 * length(model@diffy))
  } else
    pg <- model@ssd

  if (sgf[1] == 0)
  {
    sgf <- sgf[-1]
    sgf.d1 <- sgf.d1[-1,]
    sgf.d2 <- sgf.d2[,,-1]
    pg <- pg[-1]
  }

  if (hessian || infomat || modcovgrad)
    sgf.sq <- sgf^2
  if (gradient || hessian)
  {
    pgog <- as.vector(pg/sgf)
    pipgog <- pi * pgog
    if (gradient)
      gd1og <- as.matrix(sgf.d1 / sgf)
  }

  # first order derivatives

  if (gradient || (!is.null(model@transPars) && any(hessian, infomat, modcovgrad)))
  {
    d10 <- -0.5 * colSums((2 * pipgog - 1) * gd1og)

    if (!is.null(model@transPars) && version == "2") {
      d1 <- d10 * dtrans$gradient
    } else d1 <- d10

    if (barrier$mu != 0)
    {      
      if (!is.null(model@lower)) {
        d1 <- d1 + barrier$mu * bar$dl1
      }
      if (!is.null(model@upper)) {
        d1 <- d1 + barrier$mu * bar$du1
      }
    }

    # xreg

    if (!is.null(model@xreg))
    {
      if (is.list(xreg))
      {
        if (is.null(xreg$dxreg)) {
          dxreg <- model@fdiff(model@xreg, frequency(model@y))
        } else
          dxreg <- xreg$dxreg

        if (is.null(xreg$fft.dxreg)) {
          fft.dxreg <- apply(dxreg, 2, fft)
        } else 
          fft.dxreg <- xreg$fft.dxreg
      } else {
        dxreg <- model@fdiff(model@xreg, frequency(model@y))
        fft.dxreg <- apply(dxreg, 2, fft)
      }

      tmp <- as.vector(fft(model@diffy) - fft.dxreg %*% cbind(xregcoefs))
      #fft(diff(model@y - xcoef * xreg, 4))
      # "sgf" is independent of xreg
      part <- -Re(fft.dxreg) * Re(tmp) - Im(fft.dxreg) * Im(tmp)
      #part <- part / length(model@diffy) #(pi2 * length(model@diffy))
      #d1xreg <- 2 * pi * sum(part / sgf)
      d1xreg <- colSums(part / sgf) / n
      names(d1xreg) <- xregnms
      d1 <- c(d1, d1xreg)
    }

  } else d1 <- NULL

  # second order derivatives

  if (hessian)
  {
    tmp1 <- (0.5 - 2 * pipgog) / sgf.sq
    gd1cp <- apply(sgf.d1, MARGIN = 1, FUN = tcrossprod)
    d2 <- matrix(-colSums(tmp1 * t(gd1cp)), nrow = length(model@pars) - ncxreg)
    tmp2 <- (-0.5 + pipgog) / sgf

if (model@model %in% c("cycle", "trend-cycle")) {
  tmp <- matrix(0, nrow(d2), nrow(d2))
  for(i in seq(dim(sgf.d2)[3]))
    tmp <- tmp + tmp2[i] * sgf.d2[,,i]
  d2 <- d2 - tmp
} else {
  ##NOTE
  # sgf.d2 contains zeros in pure variance model
  # check it again before removing this line (tmp2 obtained above is not needed)
  diag(d2) <- diag(d2) - colSums(tmp2 * sgf.d2)
}

    if (!is.null(model@transPars) && version == "2")
      d2 <- tcrossprod(dtrans$gradient) * d2 + d10 * dtrans$hessian

    if (!is.null(model@xreg))
    {
      d2varxreg <- NULL
      for (i in seq_len(ncxreg))
      {
        tmp <- -colSums(sgf.d1 * part[,i] / sgf.sq) / n
        d2varxreg <- rbind(d2varxreg, tmp)
      }

      d2xregdiag <- colSums(Mod(fft.dxreg)^2 / sgf) / n

      if (ncxreg > 1)
      {
        d2xregxreg <- rep(0, ncxreg - 1)
        for (i in seq_len(ncxreg)) for (j in seq_len(ncxreg))
        {
          if (i > j)
            d2xregxreg[i-1] <- sum((Re(fft.dxreg[,i]) * Re(fft.dxreg[,j]) + 
              Im(fft.dxreg[,i]) * Im(fft.dxreg[,j])) / sgf) / n
        }

        d2x.aux <- diag(d2xregdiag)
        d2x.aux <- diag(d2xregdiag)
        d2x.aux[lower.tri(d2x.aux)] <- d2xregxreg
        d2x.aux[upper.tri(d2x.aux)] <- t(d2x.aux)[upper.tri(d2x.aux)]
        d2x.aux <- rbind(t(d2varxreg), d2x.aux)
      } else {
        d2x.aux <- c(d2varxreg, d2xregdiag)
      }

      d2 <- cbind(rbind(d2, d2varxreg), d2x.aux)
    }

    rownames(d2) <- colnames(d2) <- names(model@pars)

    ##NOTE see if any changes would be required if xreg is not NULL
    if (barrier$mu != 0)
    {
      if (!is.null(model@lower))
        diag(d2) <- diag(d2) + barrier$mu * bar$dl2
      if (!is.null(model@upper))
        diag(d2) <- diag(d2) + barrier$mu * bar$du2
    }

  } else d2 <- NULL

  if (infomat || modcovgrad)
  {
    if (hessian)
    {
      # ncol(sgf.d1) or length(model@pars) - ncxreg
      gcov <- matrix(0.5 * colSums(t(gd1cp) / sgf.sq), nrow = ncol(sgf.d1))
    } else
    {
      gcov <- 0
      for (i in seq(along = sgf))
        gcov <- gcov + tcrossprod(sgf.d1[i,]) / sgf.sq[i]
      gcov <- 0.5 * gcov
    }

    if (!is.null(model@transPars) && version == "2") {
      if (modcovgrad) 
        gcovmod <- gcov
      gcov <- tcrossprod(dtrans$gradient) * gcov
    } else if (modcovgrad) gcovmod <- gcov

    if (barrier$mu != 0)
    {
      if (!is.null(model@lower))
        diag(gcov) <- diag(gcov) + barrier$mu * bar$dl2
      if (!is.null(model@upper))
        diag(gcov) <- diag(gcov) + barrier$mu * bar$du2
    }

    if (!is.null(model@xreg))
    {
      ##NOTE
      # the Hessian is used, the information matrix for the 
      # coefficients of regressor would depend on the spectral 
      # generating function of the regressors 

      if (!hessian)
      {
        d2varxreg <- NULL
        for (i in seq_len(ncxreg))
        {
          tmp <- -colSums(sgf.d1 * part[,i] / sgf.sq) / n
          d2varxreg <- rbind(d2varxreg, tmp)
        }

        d2xregdiag <- colSums(Mod(fft.dxreg)^2 / sgf) / n

        if (ncxreg > 1)
        {
          d2xregxreg <- rep(0, ncxreg - 1)
          for (i in seq_len(ncxreg)) for (j in seq_len(ncxreg))
          {
            if (i > j)
              d2xregxreg[i-1] <- sum((Re(fft.dxreg[,i]) * Re(fft.dxreg[,j]) + 
                Im(fft.dxreg[,i]) * Im(fft.dxreg[,j])) / sgf) / n
          }
          d2x.aux <- diag(d2xregdiag)
          d2x.aux <- diag(d2xregdiag)
          d2x.aux[lower.tri(d2x.aux)] <- d2xregxreg
          d2x.aux[upper.tri(d2x.aux)] <- t(d2x.aux)[upper.tri(d2x.aux)]
          d2x.aux <- rbind(t(d2varxreg), d2x.aux)
        } else {
          d2x.aux <- c(d2varxreg, d2xregdiag)
        }
      }

      gcov <- cbind(rbind(gcov, d2varxreg), d2x.aux)
      rownames(gcov) <- colnames(gcov) <- names(model@pars)
    }

  } else gcov <- NULL

  if (modcovgrad)
  {
    if (!is.null(model@xreg))
      stop(sQuote("modcovgrad"), 
      " is not available with non-null regressors ", sQuote("model@xreg"))

    if (!hessian)
      tmp2 <- (-0.5 + pipgog) / sgf

    if (model@model %in% c("cycle", "trend-cycle"))
    { # for devel version
      tmp <- matrix(0, nrow(d2), nrow(d2))
      for(i in seq(dim(sgf.d2)[3]))
        tmp <- tmp + tmp2[i] * sgf.d2[,,i]
      gcovmod <- gcovmod + tmp
    } else {
        diag(gcovmod) <- diag(gcovmod) + colSums(tmp2 * sgf.d2)
    }

    if (!is.null(model@transPars) && version == "2")
      gcovmod <- tcrossprod(dtrans$gradient) * gcovmod - d10 * dtrans$hessian

    if (barrier$mu != 0)
    {
      if (!is.null(model@lower))
        diag(gcovmod) <- diag(gcovmod) + barrier$mu * bar$dl2
      if (!is.null(model@upper))
        diag(gcovmod) <- diag(gcovmod) + barrier$mu * bar$du2
    }
  } else gcovmod <- NULL

  list(gradient = d1, hessian = d2, infomat = gcov, modcovgrad = gcovmod)
}

mloglik.fd.grad <- function(x, model, xreg = NULL,
  barrier = list(type = c("1", "2"), mu = 0),
  inf)
#argument 'inf' is not used here but it is needed when used within 'optim'
#'maxlik.fd.optim' where this function is passed as the gradient
{
  if (!missing(x)) if (!is.null(x))
    model@pars <- x

  mloglik.fd.deriv(model = model, xreg = xreg, 
    gradient = TRUE, hessian = FALSE, infomat = FALSE, modcovgrad = FALSE,
    barrier = barrier, version = "2")$gradient
}
