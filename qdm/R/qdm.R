# qdm.R
#
# last mod: Aug/21/2014, NU


## QDM function
qdmfun <- function(x, y, p, response = c("logistic", "guessing", "gumbel",
                   "gompertz", "weibull", "cauchy", "shepardA", "shepardAneg",
                   "shepardB", "shepardBneg", "shepardD", "shepardDneg",
                   "shepardE", "shepardEneg", "shepardF", "shepardFneg"),
                   bias = 0) {

  s <- p[1] + p[2]*x + p[3]*x^2 + 2*p[4]*abs((x-bias) - y) + p[5]*y + p[6]*y^2

  response <- match.arg(response)
  switch(EXPR = response,
     logistic = 1/(1 + exp(-p[7]*s - p[8])),
     guessing = p[9] + (1 - p[9])*1/(1 + exp(-p[7]*s - p[8])),
       gumbel = exp(-exp((-p[7]*s - p[8])/p[9])),
     gompertz = p[7]*exp(-p[8]*exp(-p[9]*s)),
      weibull = pweibull(s, shape=p[7], scale=p[8]),
       cauchy = pcauchy(s, location=p[7], scale=p[8]),
     shepardA = 1 - s*p[7] + (s*p[7])*log(abs(s*p[7])),
  shepardAneg = 1 - (1 - s*p[7] + (s*p[7])*log(abs(s*p[7]))),
     shepardB = 1 - (s*p[7])^2 + p[8]*(s*p[7])*log(abs(s*p[7])),
  shepardBneg = 1 - (1 - (s*p[7])^2 + p[8]*(s*p[7])*log(abs(s*p[7]))),
     shepardD = 1 - p[8]*(s*p[7]) + (s*p[7])^2,
  shepardDneg = 1 - (1 - p[8]*(s*p[7]) + (s*p[7])^2),
     shepardE = 1 - p[8]*s*p[7] + p[8]*(s*p[7])^2 - (s*p[7])^3,
  shepardEneg = 1 - (1 - p[8]*s*p[7] + p[8]*(s*p[7])^2 - (s*p[7])^3),
     shepardF = exp(-p[7]*s - p[8]),
  shepardFneg = 1 - exp(-p[7]*s - p[8])
  )
}

## QDM objective function
objfun <- function(p, psi, estimfun = c("minchi2", "ols", "wls"), ...){
  yhat <- outer(psi$x, psi$y, function(x, y) qdmfun(x, y, p, ...))

  estimfun <- match.arg(estimfun)
  out <- switch(EXPR = estimfun,
    minchi2 =
      sum(na.omit(as.vector((psi$freq - psi$ntrials*yhat)^2 /
                            (psi$ntrials*yhat * (1 - yhat))))),
    ols =
      sum((psi$prob - yhat)^2),
    wls = {
      weights <- 1 / (psi$prob*(1 - (psi$prob - 0.01)))
      sum(na.omit(as.vector(weights*((psi$freq - psi$ntrials*yhat)^2))))
    }
  )
  #yhat[is.na(yhat)] <- 0.5
  if (any(yhat < 0) | any(yhat > 1)) Inf else out
}


## QDM user interface
qdm <- function(psi, start, respfun = c("logistic", "guessing", "gumbel",
                   "gompertz", "weibull", "cauchy", "shepardA", "shepardAneg",
                   "shepardB", "shepardBneg", "shepardD", "shepardDneg",
                   "shepardE", "shepardEneg", "shepardF", "shepardFneg"),
                   bias = 0, estimfun = c("minchi2", "ols", "wls"),
                   optimizer = c("optim", "nlm"), optimargs = list()){

  stopifnot(class(psi) == "psi")
  respfun <- match.arg(respfun)
  estimfun <- match.arg(estimfun)

  ## Optimize
  optimizer <- match.arg(optimizer)
  if (optimizer == "nlm") {
    optArgs <- list(f=objfun, p=start, psi=psi, estimfun=estimfun,
                    response=respfun, bias=bias)
    optArgs <- c(optArgs, as.list(optimargs))
    optimout <- do.call(nlm, optArgs)
    coefficients <- optimout$estimate
  } else {  # optim()
    optArgs <- list(par=start, fn=objfun, psi=psi, estimfun=estimfun,
                    response=respfun, bias=bias)
    optArgs <- c(optArgs, as.list(optimargs))
    optimout <- do.call(optim, optArgs)
    coefficients <- optimout$par
  }

  retval = list(optimout=optimout, coefficients=coefficients, psi=psi,
                respfun=respfun, bias=bias)
  class(retval) <- "qdm"
  retval
}
# TODO Give warning if optimizer does not give new parameters?


print.qdm <- function(x, digits = max(3L, getOption("digits") - 3L),
                      ...){
  cat("\nQuadrilateral dissimilarity models\n\n")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L, 
                quote = FALSE)
  cat("\n")
  invisible(x)
}


## FIX ME: args to qdmfun
predict.qdm <- function(object, x = object$psi$x, y = object$psi$y,
                        respfun = object$respfun, bias = object$bias, ...){
  yhat <- outer(x, y, function(x, y) qdmfun(x, y, coef(object), respfun, bias))
  dimnames(yhat) <- list(x, y)
  yhat
}

persp.qdm <- function(x, col="gray", zlim=0:1, phi=10, theta=-25,
                      xlab="OA1", ylab="OA2", zlab="Predictions", ...){
  persp(x$psi$x, x$psi$y, predict(x), col=col, zlim=zlim, phi=phi,
        theta=theta, xlab=xlab, ylab=ylab, zlab=zlab, ...)
  }

