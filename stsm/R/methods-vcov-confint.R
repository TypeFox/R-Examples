
#NOTE do not pass 'vcov.type[1]' or 'type[1]' since method 'vcov()' 
#uses length(vcov.type) to see whether a choice was made or not and 
#set default value for this argument

##FIXME update with xreg, see mloglik.td.deriv

vcov.stsmFit <- function(object, 
  type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), ...)
  #domain = c("frequency", "time"), 
  #td.args = list(), ...)
{ 
  #domain <- match.arg(domain)[1]
  domain <- ifelse (grepl("maxlik.fd", object$call[1]), 
    "frequency", "time")

  # if 'length(type) > 1' no choice is made, the default is taken
  if (length(type) > 1)
  { 
    if (domain == "time") { # default for the time domain
      type <- "infomat"
    } else # default for the frequency domain
      type <- "hessian"
  } else 
    type <- match.arg(type)[1]

  if (type == "optimHessian")
  {
    if (!is.null(object$hessian)) {
      res <- solve(object$hessian)
      nms <- c(names(object$model@pars), names(object$xreg$coef))
      rownames(res) <- colnames(res) <- nms
      return(res)
    } else {
      type <- "hessian"
      warning("The model was not fit by ", sQuote("optim"), 
      #"or ", sQuote("hessian"), " was set to FALSE",
      ", hence, ", sQuote("type = 'optimHessian'"), " cannot be applied. ",
      "Argument ", sQuote("type"), " was changed to ", sQuote("hessian"), ".")
    }
  } 
#   else if (type == "infomat")
#   {
#     if (!is.null(object$infomat)) {
#       res <- solve(object$hessian)
#       nms <- c(names(object$model@pars), names(object$xreg$coef))
#       rownames(res) <- colnames(res) <- nms
#       return(res)
#     }
#   }

  dogr <- type %in% c("OPG", "sandwich")
  dohe <- type %in% c("hessian", "sandwich")
  doim <- type == "infomat"

  switch(domain,
    "frequency" = {
      if (is.null(object$model@cpar)) {
        res <- mloglik.fd.deriv(model = object$model, gradient = dogr, 
          hessian = dohe, infomat = doim, modcovgrad = FALSE,
          version = "2")
      } else 
        res <- mcloglik.fd.deriv(model = object$model, 
          gradient = dogr, hessian = dohe, infomat = doim)
    },

    "time" = {
      if (type == "hessian" || type == "sandwich")
      {
        type <- "infomat"
        doim <- TRUE #type <- "infomat"
        warning("The analytical Hessian is not available for the time domain version. ",
        "The information matrix was used instead, changed to ", 
        sQuote("type = 'infomat'"), ".")
      }

      ##NOTE 
      # P0cov = FALSE is considered regardless of the value 
      # in the optimization procedure
      P0cov <- FALSE

      res <- mloglik.td.deriv(model = object$model, 
        gradient = dogr, infomat = doim,
        KF.args = list(P0cov = P0cov), version = "1", kfres = NULL)
    })

  switch(type, 
    "hessian" = res <- solve(res$hessian), 
    "infomat" = res <- solve(res$infomat),
    "OPG" = res <- solve(crossprod(rbind(res$gradient))), 
    "sandwich" = {
      invHess <- solve(res$hessian)
      opg <- crossprod(rbind(res$gradient))
      res <- invHess %*% opg %*% invHess
    })
  rownames(res) <- colnames(res) <- names(object$model@pars)

  res
}

vcov.stsm <- function(object, 
  type = c("hessian", "infomat", "OPG", "sandwich"), #"optimHessian"
  domain = c("frequency", "time"), ...)
{
  #type = "optimHessian" does not apply here

  # artificial 'stsmFit' object with the elements requied by 'vcov.stsmFit'
  mcall <- switch(match.arg(domain)[1], 
    "frequency" = "maxlik.fd", time = "maxlik.td")

  x <- list(call = mcall, model = object)
  class(x) <- "stsmFit"

  vcov.stsmFit(x, type = type, ...)
}

confintStsmFitVcov <- function(object, parm, level = 0.95, 
  type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), ...)
{
# adapted from 'confint' in 'stats' package

  # from 'stats' package file 'confint.R' (added 'sep = ""')
  format.perc <- function(probs, digits)
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, 
      digits = digits), "%", sep = "")

  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) parm <- pnames
  else if(is.numeric(parm)) parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object, type = type, ...)))[parm]
  ci[] <- cf[parm] + ses %o% fac

  # enforce bounds
  lowerb <- object$model@lower[parm]
  upperb <- object$model@upper[parm]
  if (any(ci < lowerb) || any(ci > upperb))
  {
    ci <- pmax(ci, lowerb)
    ci <- pmin(ci, upperb)  
    warning(sQuote("confint"), " did not account for the fact that ", 
    "some of the parameters lied\n outside the boundary of the parameter space.")
  }

  ci
}

confintStsmFitBoot <- function(object, parm, level = 0.95, breps = 100, ...)
{
##NOTE if length(parm)==1 the other parameters could be fixed to the
#values from object@model (local optimum) in the bootstrapped model 'mb'

  pnames <- names(coef(object))
  if (missing(parm)) parm <- pnames
  else if(is.numeric(parm)) parm <- pnames[parm]

  ref <- ifelse(length(object$model@diffy) %% 2 == 0, 1, 0)

  nh <- floor(length(object$model@diffy) / 2)
  nhmref <- nh - ref

  nhp1 <- nh + 1
  #nhp2 <- nh + 2
  pi4 <- 4 * pi
  npars <- length(object$model@pars)
  bpars <- matrix(nrow = breps, ncol = npars)
  colnames(bpars) <- parm
  mb <- object$model
  #mb@y <- ts(seq_along(object$model@y))
  #attributes(mb@y) <- attributes(object$model@y)
  #mb@diffy <- ts(seq_along(object$model@diffy))

  sgf <- stsm.sgf(object$model, FALSE, FALSE, FALSE)
  if (is.null(mb@sgfc))
    mb@sgfc <- sgf$constants
  sgf <- sgf$sgf[seq_len(nhp1)]

  # bootstrapped periodogram based on Dahlhaus-Janas:96

  for (j in seq_len(breps))
  {
    x <- c(rchisq(1, df = 1), rchisq(nhmref, df = 2))
    if (ref == 1)
      x <- c(x, rchisq(1, df = 1))
    pgb <- x * sgf
    pgb <- sgf * x / pi4
    pgb[1] <- pgb[1] * 2
    if (ref == 1)
      pgb[nh] <- pgb[nh] * 2
    if (ref == 1) {
      mb@ssd <- c(pgb, rev(pgb[-c(1, nhp1)]))
    } else 
      mb@ssd <- c(pgb, rev(pgb[-1]))

    resb <- maxlik.fd.scoring(m = mb, step = NULL, information = "expected", 
      ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
      barrier = list(mu = 0), debug = FALSE,
      control = list(maxit = 100, tol = 0.001, trace = TRUE))

    bpars[j,] <- coef(resb)
  }

  a <- (1 - level)/2
  a <- c(a, 1 - a)
  ci <- t(apply(data.matrix(bpars[,parm]), MARGIN = 2, FUN = quantile, probs = a))
  rownames(ci) <- parm

  ci
}

confint.stsmFit <- function(object, parm, level = 0.95, 
  type = c("vcov", "bootstrap"),
  vcov.type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), 
  breps = 100, ...)
{
  type <- match.arg(type)[1]
  switch(type, 
    "vcov" = confintStsmFitVcov(object, parm, level, vcov.type),
    "bootstrap" = confintStsmFitBoot(object, parm, level, breps))
}
