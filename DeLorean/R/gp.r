#' Centralises a periodic position into [period/2, period)
#' by shifting by n*period, where n is an integer
#'
#' @param x Position
#' @param period Period
#'
centralise <- function(x, period=1) {
  y <- x / period + .5
  (y - floor(y) - .5) * period
}

#' Makes a distance periodic
#'
#' @param r Distance
#' @param period The period
#'
cov.periodise <- function(r, period) {
    period * sin(r * pi / period) / 2;
}

#' Matern 3/2 covariance function
#'
#' @param r Distance
#' @param l Length scale
#'
cov.matern.32 <- function(r, l) {
    x <- sqrt(3) * abs(r / l);
    (1 + x) * exp(- x);
}

#' Calculate distances between vectors of time points
#'
#' @param tau.1 First vector of time points
#' @param tau.2 Second vector of time points (defaults to first if not given)
#' @param period Period if periodic
#'
cov.calc.dists <- function(tau.1, tau.2=tau.1, period=NULL) {
    # Calculate the distances
    r <- outer(tau.1, tau.2, FUN="-")
    # Make periodic if necessary
    if (! is.null(period) && period > 0) {
        r <- cov.periodise(r, period)
    }
    r
}


#' Calculate distances over estimated pseudotimes and
#' test inputs.
#'
#' @param dl de.lorean object
#' @param tau The pseudotimes to use
#' @param include.test Also include the pseudotimes for the test inputs
#'
cov.calc.dl.dists <- function(dl,
                              tau=tau.for.sample(dl),
                              include.test=TRUE) {
    with(dl, {
        if (length(tau) == 1 && "capture" == tau) {
            tau <- dl$cell.map$obstime
        }
        if (include.test) {
            tau <- c(tau, test.input)
        }
        if (exists('periodic', opts) && opts$periodic) {
            cov.calc.dists(tau, period=opts$period)
        } else {
            cov.calc.dists(tau)
        }
    })
}

#' Calculate covariance structure for gene over pseudotimes and test
#' inputs.
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param tau The pseudotimes to use
#' @param include.test Also include the pseudotimes for the test inputs
#' @param psi Temporal variation
#' @param omega Noise
#'
cov.calc.gene <- function(dl,
                          gene.idx,
                          cov.fn=cov.matern.32,
                          tau=tau.for.sample(dl),
                          include.test=TRUE,
                          psi = sampled.gene.param(dl, gene.idx, "psi"),
                          omega=sampled.gene.param(dl, gene.idx, "omega"))
{
    with(dl, {
        r <- cov.calc.dl.dists(dl,
                               tau=tau,
                               include.test=include.test)
        (
            psi * cov.fn(r, opts$length.scale)  # Structure
            + omega * diag(nrow(r))  # Noise
        )
    })
}

#' Calculate covariance for gene over test inputs when conditioned on
#' data at estimated pseudotimes.
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param tau The pseudotimes to use
#'
cov.calc.gene.conditioned <- function(dl,
                                      gene.idx,
                                      cov.fn=NULL,
                                      tau=tau.for.sample(dl))
{
    if (is.null(cov.fn)) {
        cov.fn <- cov.matern.32
    }
    with(dl, {
        Sigma <- cov.calc.gene(dl, gene.idx, cov.fn, tau=tau)
        num.cells <- nrow(dl$cell.map)
        # From Appendix A.2 of Rasmussen/Williams GP book
        slice.obsd <- 1:num.cells
        slice.test <- (num.cells+1):nrow(Sigma)
        .B <- Sigma[slice.obsd,slice.obsd]
        .A <- Sigma[slice.test,slice.test]
        .C <- Sigma[slice.test,slice.obsd]
        .A - .C %*% qr.solve(.B, t(.C))
    })
}

#' Calculate covariances for all genes when conditioned on data at
#' estimated pseudotimes.
#'
#' @param dl de.lorean object
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param tau The pseudotimes to use
#'
cov.all.genes.conditioned <- function(dl,
                                      cov.fn=NULL,
                                      tau=tau.for.sample(dl))
{
    vapply(1:nrow(dl$gene.map),
           function(gene.idx) cov.calc.gene.conditioned(dl,
                                                        gene.idx,
                                                        cov.fn,
                                                        tau=tau),
           FUN.VALUE=matrix(0,
                            nrow=length(dl$test.input),
                            ncol=length(dl$test.input)))
}

#' Add posterior representation to a plot.
#'
#' @param gp Plot object
#' @param .data Data frame containing variables to plot (mean, var)
#'        phi, predictedvar)
#' @param color Color to use
#' @param line.alpha Alpha to use for mean line
#' @param ribbon.alpha Alpha to use for variance ribbon
#'
plot.add.mean.and.variance <- function(gp,
                                       .data=NULL,
                                       color='black',
                                       line.alpha=.3,
                                       ribbon.alpha=.1) {
    (gp
        + geom_line(data=.data,
                    aes(x=x, y=mean),
                    color=color,
                    alpha=line.alpha)
        + geom_ribbon(data=.data,
                      aes(x=x,
                          ymin=mean-2*sqrt(var),
                          ymax=mean+2*sqrt(var)),
                      fill=color,
                      alpha=ribbon.alpha))
}

#' The log marginal likelihood. See "2.3 Varying the Hyperparameters"
#' on page 19 of Rasmumssen and Williams' book.
#'
#' @param y The targets.
#' @param K The covariance matrix (kernel), not needed if U is provided.
#' @param U Cholesky decomposition of K (chol(K)).
#'
gp.log.marg.like <- function(y, K=NULL, U=chol(K)) {
    alpha <- backsolve(U, backsolve(U, y, transpose = TRUE))
    -(
        (t(y) %*% alpha) / 2
        + sum(diag(U))
        + length(y) * log(2 * pi)
    )
}

#' Predictive mean, variance and log marginal likelihood of a GP.
#' See "2.3 Varying the Hyperparameters"
#' on page 19 of Rasmumssen and Williams' book.
#'
#' @param y The targets.
#' @param K The covariance matrix (kernel) for input points, not needed if U is provided.
#' @param Kstar The cross covariance matrix (kernel)
#' @param Kstarstar The cross covariance matrix (kernel) for test points
#' @param U Cholesky decomposition of K
#'
gp.predict <- function(y, K=NULL, Kstar, Kstarstar, U=chol(K)) {
    alpha <- backsolve(U, backsolve(U, y, transpose = TRUE))
    v <- backsolve(U, Kstar, transpose = TRUE)
    list(
        mu = t(Kstar) %*% alpha,  # fstar
        Sigma = Kstarstar - t(v) %*% v,  # Vstar
        logpy = -(
            (t(y) %*% alpha) / 2
            + sum(diag(U))
            + length(y) * log(2 * pi)
        )
    )
}

#' Convert the output of gp.predict() into a data.frame.
#'
#' @param predictions The predictions
#'
gp.predictions.df <- function(predictions) {
    with(predictions, data.frame(mu=mu,
                                 Sigma=diag(Sigma),
                                 idx=1:length(mu)))
}

#' Condition a Guassian on another.
#' See Eqn. A.6
#' on page 200 of Rasmumssen and Williams' book.
#'
#' @param y y
#' @param mu.x Mean of x
#' @param mu.y Mean of y
#' @param .A Var(X)
#' @param .B Var(Y)
#' @param .C Cov(X, Y)
#' @param U Cholesky decomposition of .B
#'
gaussian.condition <- function(
    y,
    .A,
    .B,
    .C,
    mu.x=rep(0, nrow(.A)),
    mu.y=rep(0, nrow(.B)),
    U=chol(.B))
{
    alpha <- backsolve(U, t(.C))
    print(dim(alpha))
    print(dim(.C))
    print(dim(alpha %*% .C))
    print(dim(.A))
    list(mu = mu.x + as.vector((y - mu.y) %*% alpha),
         Sigma = .A - .C %*% alpha)
}

#' The expected within sample variance of a Gaussian with the given covariance.
#'
#' @param K Covariance
#'
#' @export
#'
expected.sample.var <- function(K) mean(diag(K)) - mean(K)

#' Calculate the covariance structure of the tau
#'
#' @param cov.fn Covariance function
#' @param length.scale Length scale
#' @param period The period of the covariance function
#' @param tau Pseudotimes
#' @param psi The temporal variation for each gene
#' @param omega The noise level for each gene
#'
#' @export
#'
gene.covariances <- function(
  cov.fn,
  length.scale,
  period=NULL,
  tau,
  psi,
  omega)
{
  stopifnot(length(psi) == length(omega))
  #
  # Calculate covariance between inducing inputs
  r <- cov.calc.dists(tau.1=tau, period=period)
  within(list(), {
    K <- cov.fn(r, length.scale)
    .G <- length(psi)
    .T <- length(tau)
    gene.K <- array(dim=c(.G, .T, .T))
    gene.K.chol <- array(dim=c(.G, .T, .T))
    log.K.det <- rep(NA, .G)
    for (g in 1:.G) {
      gene.K[g,,] <- psi[g] * K + diag(omega, .T)
      gene.K.chol[g,,] <- chol(gene.K[g,,])
      log.K.det[g] <- sum(log(diag(gene.K.chol[g,,])))
    }
  })
}

#' Calculate the covariance structure of evenly spread tau and
#' create a function that calculates the log likelihood of
#' orderings.
#'
#' @param dl The DeLorean object
#' @param cov.fn The covariance function
#'
#' @export
create.ordering.ll.fn <- function(dl, cov.fn=cov.fn.for(dl)) within(dl, {
  if (! exists('even.tau.cov')) {
    even.tau.cov <- gene.covariances(cov.fn, opts$length.scale,
                                    opts$period, even.tau.spread(dl),
                                    gene.var$psi.hat,
                                    gene.var$omega.hat)
    expr.centre <- t(scale(t(expr), scale=FALSE, center=TRUE))
    ordering.ll <- function(o) {
      # Return the sum of the marginal log likelihoods for each gene
      sum(sapply(1:stan.data$G,
                 function(g) {
                   gp.log.marg.like(expr.centre[g,o],
                                    U=even.tau.cov$gene.K.chol[g,,])
                 }))
    }
    order.inits <- list()
  }
})

#' Calculate the covariance structure of the inducing points
#'
#' @param cov.fn Covariance function
#' @param length.scale Length scale
#' @param period The period of the covariance function
#' @param tau Pseudotimes
#' @param num.inducing How many inducing points to use
#' @param u The inducing points
#'
#' @export
#'
inducing.covariance <- function(
  cov.fn,
  length.scale,
  period=NULL,
  tau,
  num.inducing,
  u=seq(min(tau), max(tau), length.out=num.inducing))
{
  #
  # Calculate covariance between inducing inputs
  r.uu <- cov.calc.dists(tau.1=u, period=period)
  K.uu <- cov.fn(r.uu, length.scale)
  within(list(), {
    u <- u
    K.uu.chol <- chol(K.uu)
    .T <- length(tau)
    .U <- length(u)
    #
    # Calculate approximate covariance between evenly spread pseudotimes
    K.ut <- cov.calc.dists(tau.1=u, tau.2=tau, period=period)
    K.tau.sqrt <- backsolve(K.uu.chol, K.ut, transpose=TRUE)
  })
}

# Calculate gene specific approximate precision from inducing covariance
# structure
#
induced.gene.prec <- function(inducing, psi, omega) with(inducing, {
  # print(.U)
  # print(dim(inducing$K.tau.sqrt))
  # print(dim(diag(omega**2/psi, .U)))
  # print(dim(omega*tcrossprod(K.tau.sqrt)))
  chol.term <- chol(diag(omega**2/psi, .U) + omega*tcrossprod(K.tau.sqrt))
  print('Done chol')
  # print(dim(crossprod(backsolve(chol.term, K.tau.sqrt, transpose=TRUE))))
  K.inv <-
    diag(1/omega, .T) -
    crossprod(backsolve(t(chol.term), K.tau.sqrt, transpose=TRUE))
})
