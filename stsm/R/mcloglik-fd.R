
##FIXME TODO add "xreg" functionality (currently only in the prototype)

mcloglik.fd <- function(x, model, xreg = NULL,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999)
{
  if (!missing(x)) if (!is.null(x))
    model@pars[] <- x

  n <- length(model@diffy)
  nh <- n/2
  pi2 <- 2 * pi

  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (pi2 * n)
  } else
    pg <- model@ssd

  sgf <- stsm.sgf(model, FALSE, FALSE, FALSE)$sgf
  cvar <- (pi2 / n) * sum(pg / sgf)

  mll <- nh * (log(pi2) + 1) + nh * log(cvar) + 0.5 * sum(log(sgf))

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

mcloglik.fd.deriv <- function(model, xreg = NULL, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE)
{
##TODO 'dtpars'
  if (!is.null(model@transPars))
    stop(sQuote("mcloglik.fd.deriv"), "is not implemented for model with non-null", 
      sQuote("model@transPars"), ".")

  g <- h <- im <- NULL
  npars <- length(model@pars)
  n <- length(model@diffy)
  nh <- n / 2
  pi2 <- 2 * pi

  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (pi2 * n)
  } else
    pg <- model@ssd

  tmp <- stsm.sgf(model, FALSE, FALSE, FALSE)
  part <- data.matrix(tmp$constants)
  sgf <- tmp$sgf
  sgf.sq <- sgf^2

  ref <- match(names(model@cpar), colnames(part))
  part <- data.matrix(part[,-ref])

  # gradient

  if (gradient)
  {
    dvc <- colSums(part * (c(pg) / sgf.sq))
    g <- -nh * dvc / sum(pg / sgf) + 0.5 * colSums(part / sgf)
  }

  # hessian 

  if (hessian)
  {
    h <- matrix(nrow = npars, ncol = npars)
    for (i in seq(npars))
    {
      for (j in seq(npars))
      {
        num <- sum(part[,i] * pg / sgf.sq)
        dem <- sum(pg / sgf)
        aux <- -2 * (part[,i] * part[,j] * pg / sgf^3) * dem + num * part[,j] * pg / sgf.sq
        h[i,j] <- -nh * sum(aux) / dem^2 - 0.5 * sum(part[,i] * part[,j] /sgf.sq)
      }
    }
  }

  # IM

  if (infomat)
  {
    im <- matrix(nrow = npars, ncol = npars)
    ipi2sgf <- 1 / (pi2 * sgf)
    ipi2sgfsq <- 1 / (pi2 * sgf.sq)
    ni2pi <- n / pi2

    for (i in seq(npars))
    {
      for (j in seq(npars))
      {
        aux1 <- - 2 * sum(part[,i] * part[,j] * ipi2sgfsq) * ni2pi
        aux2 <- sum(part[,i] * ipi2sgf) * sum(part[,j] * ipi2sgf)
        aux3 <- 0.5 * sum(part[,i] * part[,j] / sgf.sq)
        im[i,j] <- -(nh * ((aux1 + aux2) / ni2pi^2) + aux3)
      }
    }
  }

  list(gradient = g, hessian = h, infomat = im)
}

mcloglik.fd.grad <- function(x, model, xreg = NULL, inf, barrier)
#arguments 'inf' and 'barrier' are not used here but it is needed when used within optim 
#'maxlik.fd.optim' where this function is passed as the gradient
{
  if (!missing(x)) if (!is.null(x))
    model@pars <- x

  mcloglik.fd.deriv(model = model, 
    gradient = TRUE, hessian = FALSE, infomat = FALSE)$gradient
}
