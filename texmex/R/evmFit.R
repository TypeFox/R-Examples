evmFit <- function(data, family, ...,
                   prior="none", start=NULL,
                   priorParameters = NULL,
                   maxit = 10000, trace = 0, hessian=TRUE) {

  penFactory <- switch(prior,
                       laplace=.make.lasso.penalty,
                       lasso=.make.lasso.penalty,
                       l1=.make.lasso.penalty,
                       quadratic=.make.quadratic.penalty,
                       gaussian=.make.quadratic.penalty,
                       none=.make.dummy.penalty,
                       function() {stop("Bad penalty ref.")})

  prior <- penFactory(priorParameters)

  log.lik <- family$log.lik(data, ...)

  evm.lik <- function(par) {
    min(-log.lik(par), 1e6) + prior(par)
  }

  nco <- sum(sapply(data$D, ncol))

  if (is.null(start)){
    if (is.null(family$start)){ start <- runif(length(nco), -.1, .1) }
    else { start <- family$start(data) }
  }

   o <- optim(par = start, fn = evm.lik,
              control = list(maxit = maxit, trace = trace),
              hessian = hessian)

    invisible(o)
}
