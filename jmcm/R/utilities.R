#' @title Extract or Get Generalized Components from a Fitted Joint Mean
#' Covariance Model
#'
#' @description Extract (or "get") "components" - in a generalized sense - from
#' a fitted joint mean covariance model from an object of class "jmcmMod".
#'
#' @param object a fitted joint mean covariance model of class "jmcmMod", i.e.,
#' typically the result of jmcm().
#' @param name a character vector specifying the name(s) of the "component".
#'
#' When sub.num is not specified or equal to 0, possible values are:
#' \describe{
#'   \item{\code{"m"}}{a vector of number of measurement for each subject}
#'   \item{\code{"Y"}}{response vector}
#'   \item{\code{"X"}}{model matrix for mean structure}
#'   \item{\code{"Z"}}{model matrix for covariance structure (the diagonal
#'   matrix)}
#'   \item{\code{"W"}}{model matrix for covariance structure (the lower
#'   triangular matrix)}
#'   \item{\code{"theta"}}{parameter estimates of joint mean covariance model}
#'   \item{\code{"beta"}}{parameter estimates for mean structure model}
#'   \item{\code{"lambda"}}{parameter estimates for covariace structure (the
#'   diagonal matrix)}
#'   \item{\code{"gamma"}}{parameter estimates for covariance structure (the
#'   lower triangular matrix)}
#'   \item{\code{"loglik"}}{log-likelihood, except for a constant}
#'   \item{\code{"BIC"}}{Bayesian information criterion}
#'   \item{\code{"iter"}}{number of iterations until convergence}
#'   \item{\code{"triple"}}{(p, d, q)}
#' }
#'
#' When sub.num is specified, possible values are:
#' \describe{
#'   \item{\code{"m"}}{number of measurements for subject i}
#'   \item{\code{"Y"}}{response vector for subject i}
#'   \item{\code{"X"}}{model matrix of subject i for mean structure }
#'   \item{\code{"Z"}}{model matrix of subject i for covariance structure (the
#'   diagonal matrix)}
#'   \item{\code{"W"}}{model matrix of subject i for covariance structure (the
#'   lower triangular matrix)}
#'   \item{\code{"D"}}{the estimated diagonal matrix for subject i}
#'   \item{\code{"T"}}{the estimated lower triangular matrix for subject i}
#'   \item{\code{"Sigma"}}{the estimated covariance matrix for subject i}
#'   \item{\code{"mu"}}{the estimated mean for subject i}
#' }
#'
#' @param sub.num refer to i's subject
#'
#' @examples
#' fit.mcd <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
#'   triple = c(8, 1, 3), cov.method = 'mcd')
#'
#' beta <- getJMCM(fit.mcd, "beta")
#' BIC  <- getJMCM(fit.mcd, "BIC")
#' Di   <- getJMCM(fit.mcd, "D", 10)
#'
#' @export
getJMCM <- function(object, name, sub.num) UseMethod("getJMCM")

#' @describeIn getJMCM Extract or Get Generalized Components from a Fitted Joint
#' Mean Covariance Model
#' @export
getJMCM.jmcmMod <- function(object,
  name = c("m", "Y", "X", "Z", "W", "D", "T", "Sigma", "mu",
    "theta", "beta", "lambda", "gamma", "loglik", "BIC", "iter",
    "triple"),
  sub.num = 0)
{
  if(missing(name)) stop("'name' must not be missing")

  stopifnot(is(object,"jmcmMod"))

  opt     <- object@opt
  args    <- object@args
  devcomp <- object@devcomp

  if(sub.num < 0 || sub.num > length(args$m))
    stop("incorrect value for 'sub.num'")

  m = args$m
  Y = args$Y
  X = args$X
  Z = args$Z
  W = args$W
  theta  = drop(opt$par)

  if (devcomp$dims['MCD']) obj <- new("MCD", m, Y, X, Z, W)
  if (devcomp$dims['ACD']) obj <- new("ACD", m, Y, X, Z, W)
  if (devcomp$dims['HPC']) obj <- new("HPC", m, Y, X, Z, W)

  if(sub.num == 0) {
    switch(name,
      "m" = args$m,
      "Y" = args$Y,
      "X" = args$X,
      "Z" = args$Z,
      "W" = args$W,
      "theta"  = drop(opt$par),
      "beta"   = drop(opt$beta),
      "lambda" = drop(opt$lambda),
      "gamma"  = drop(opt$gamma),
      "loglik" = opt$loglik,
      "BIC"    = opt$BIC,
      "iter"   = opt$iter,
      "triple" = object$triple,
      "n2loglik" = obj$n2loglik(theta),
      "grad"     = obj$grad(theta))
  } else {
    switch(name,
      "m" = obj$get_m(sub.num),
      "Y" = obj$get_Y(sub.num),
      "X" = obj$get_X(sub.num),
      "Z" = obj$get_Z(sub.num),
      "W" = obj$get_W(sub.num),
      "D" = obj$get_D(theta, sub.num),
      "T" = obj$get_T(theta, sub.num),
      "Sigma" = obj$get_Sigma(theta, sub.num),
      "mu" = obj$get_mu(theta, sub.num))
  }
}

lagseq <- function(time)
{
  res <- NULL
  if(length(time) != 1) {
    for(i in 2:length(time)) {
      for(j in 1:(i-1))
        res <- c(res, (time[i] - time[j]))
    }
  }
  res
}

#' @title Plot Fitted Mean Curves
#'
#' @description plot fitted mean curves
#'
#' @param object a fitted joint mean covariance model of class "jmcmMod", i.e.,
#' typically the result of jmcm().
#'
#' @examples
#' cattleA <- cattle[cattle$group=='A', ]
#' fit.mcd <- jmcm(weight | id | I(ceiling(day/14 + 1)) ~ 1 | 1, data=cattleA,
#'   triple = c(8, 3, 4), cov.method = 'mcd')
#' meanplot(fit.mcd)
#'
#' @export
meanplot <- function(object)
{
  op <- par(mfrow = c(1, 1))

  opt <- object@opt
  beta <- opt$beta
  lbta <- length(beta)

  args   <- object@args
  Y <- args[["Y"]]
  time <- args[["time"]]

  ts   <- seq(min(time), max(time), length.out = 100)

  X.ts    <- NULL
  for(i in 0:(lbta-1)) X.ts    <- cbind(X.ts, ts^i)

  Yest <- drop(X.ts %*% beta)
  plot(time, Y, xlab = "Time", ylab = "Response")
  lines(ts, Yest)
}

#' @title Plot Sample Regressograms and Fitted Curves
#'
#' @description Plot the sample regressograms based on the sample covariance
#' matrix and superimpose the corresponding fitted curves to check the model
#' fitting when the longitudinal dataset is balanced.
#'
#' @param object a fitted joint mean covariance model of class "jmcmMod", i.e.,
#' typically the result of jmcm().
#' @param time a vector of obeservation time points
#'
#' @examples
#' cattleA <- cattle[cattle$group=='A', ]
#' fit.mcd <- jmcm(weight | id | I(ceiling(day/14 + 1)) ~ 1 | 1, data=cattleA,
#'   triple = c(8, 3, 4), cov.method = 'mcd')
#' regressogram(fit.mcd, time = 1:11)
#'
#' @export
regressogram <- function(object, time)
{
  debug <- 0

  op <- par(mfrow = c(1, 2))

  opt <- object@opt

  lambda <- opt$lambda
  gamma  <- opt$gamma

  llmd <- length(lambda)
  lgma <- length(gamma)

  args   <- object@args
  dims   <- object@devcomp$dims

  m <- args[["m"]]
  Y <- args[["Y"]]
  X <- args[["X"]]
  Z <- args[["Z"]]
  W <- args[["W"]]

  if (length(unique(m)) != 1)
    stop("No regressograms. Unbalanced longitudinal dataset.")

  # create a data matrix
  DataMat <- t(Y[1:m[1]])
  for(i in 2:length(m))
  {
    DataMat <- rbind(DataMat, t(Y[(sum(m[1:(i-1)])+1):sum(m[1:i])]))
  }
  dimnames(DataMat) <- NULL

  S <- cov(DataMat)  # sample covariance matrix
  R <- cor(DataMat)  # sample correlation matrix

  # FIXME: singularity check
  C <- t(chol(S))    # Cholesky factor of S
  D <- diag(diag(C))

  # transpose of matrix T in MCD
  Tt <- t(forwardsolve(C %*% diag(diag(C)^(-1)), diag(dim(D)[1])))

  # transpose of matrix L in ACD
  Lt <- t(diag(diag(C)^(-1)) %*% C)

  ts    <- seq(min(time), max(time), length.out = 100)
  tlag  <- lagseq(time)
  tslag <- seq(min(tlag), max(tlag), length.out = 100)

  Z.ts    <- NULL
  W.tslag <- NULL
  for(i in 0:(llmd-1)) Z.ts <- cbind(Z.ts, ts^i)
  for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)

  Zlmd <- Z.ts %*% lambda
  Wgma <- W.tslag %*% gamma

  # plot regressogram for MCD, ACD or HPC
  if (dims["MCD"] == 1) {
    # the first plot
    plot(time, log(diag(D)^2), xlab="Time", ylab="Log-innovat. var.")
    lines(ts, Zlmd)

    # the second plot
    phi <- -Tt[upper.tri(Tt, diag=FALSE)]
    plot(tlag, phi, xlab="Lag", ylab="Autoregres. coeffic.")
    lines(tslag, Wgma)

  } else if (dims["ACD"] == 1) {
    # the first plot
    plot(time, log(diag(D)^2), xlab="Time", ylab="Log-innovat. var.")
    lines(ts, Zlmd)

    # the second plot
    phi <- Lt[upper.tri(Lt, diag=FALSE)]
    plot(tlag, phi, xlab="Lag", ylab="MA. coeffic.")
    lines(tslag, Wgma)

  } else if (dims["HPC"] == 1) {
    # the first plot
    H <- diag(sqrt(diag(S)))
    plot(time, log(diag(H)^2), xlab="Time", ylab="Log-variance")
    lines(ts, Zlmd)

    # the second plot
    B <- t(chol(R))
    PhiMat <- matrix(0, dim(B)[1], dim(B)[2])
    for(j in 2:dim(B)[1]) {
      for(k in 1:(j-1)) {
        tmp <- 1
        if (k != 1) {
          tmp <- prod(sin(PhiMat[j, 1:(k-1)]))
        } # if
        PhiMat[j,k] <- acos(B[j, k]/tmp)
      } # for k
    } # for j
    PhiMatt <- t(PhiMat)

    phi <- PhiMatt[upper.tri(PhiMatt, diag=FALSE)]
    plot(tlag, phi, xlab="Lag", ylab="Angles")
    lines(tslag, Wgma)
  } # HPC
}

#' @title Plot Fitted Curves and Corresponding Confidence Interval using
#' bootstrapping method
#'
#' @description Plot fitted curves and corresponding 95\% confidence interval
#' using bootstrapping method.
#'
#' @param object a fitted joint mean covariance model of class "jmcmMod", i.e.,
#' typically the result of jmcm().
#' @param nboot number of the bootstrap replications.
#'
#' @examples
#' \dontrun{
#' # It may take hours for large bootstrap replications
#' fit.mcd <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data=aids,
#'   triple = c(8, 1, 3), cov.method = 'mcd', control = jmcmControl(trace=T))
#' bootcurve(fit.mcd, nboot = 1000)
#' }
#'
#' @export
bootcurve <- function(object, nboot)
{
  debug <- 0

  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

  opt <- object@opt

  theta  <- opt$par
  beta   <- opt$beta
  lambda <- opt$lambda
  gamma  <- opt$gamma

  ltht   <- length(theta)
  lbta   <- length(beta)
  llmd   <- length(lambda)
  lgma   <- length(gamma)

  args   <- object@args
  dims   <- object@devcomp$dims

  m <- args[["m"]]
  Y <- args[["Y"]]
  X <- args[["X"]]
  Z <- args[["Z"]]
  W <- args[["W"]]
  time <- args[["time"]]

  ts    <- seq(min(time), max(time), length.out = 100)
  tslag <- seq(0, max(time) - min(time), length.out = 100)

  X.ts    <- NULL
  Z.ts    <- NULL
  W.tslag <- NULL

  for(i in 0:(lbta-1)) X.ts    <- cbind(X.ts, ts^i)
  for(i in 0:(llmd-1)) Z.ts    <- cbind(Z.ts, ts^i)
  for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)

  Yest <- drop(X.ts %*% beta)
  Zlmd <- drop(Z.ts %*% lambda)
  Wgma <- drop(W.tslag %*% gamma)

  Yest.boot <- NULL
  Zlmd.boot <- NULL
  Wgma.boot <- NULL

  result <- NULL
  for(iter in 1:nboot) {
    # generate a bootstrap sample
    index <- sample(length(m), replace=T)

    # construct corresponding arguments
    m.boot <- m[index]
    Y.boot <- NULL
    X.boot <- NULL
    Z.boot <- NULL
    W.boot <- NULL
    for(i in 1:length(m)) {
      if (index[i] == 1) {
        Y.boot <- c(Y.boot, Y[1:m[1]])
        X.boot <- rbind(X.boot, X[1:m[1], ])
        Z.boot <- rbind(Z.boot, Z[1:m[1], ])

        if (m[1] != 1) {
          first <- 1
          last  <- m[1] * (m[1] - 1) / 2
          W.boot <- rbind(W.boot, W[first:last, ])
        }
      } else {
        first <- sum(m[1:(index[i]-1)]) + 1
        last  <- sum(m[1:index[i]])
        Y.boot <- c(Y.boot, Y[first:last])
        X.boot <- rbind(X.boot, X[first:last, ])
        Z.boot <- rbind(Z.boot, Z[first:last, ])

        if (m[index[i]] != 1) {
          first <- 0
          for(j in 1:(index[i]-1)) {
            first <- first + m[j] * (m[j] - 1) / 2
          }
          last  <- first  + m[index[i]] * (m[index[i]] - 1) / 2
          first <- first + 1

          W.boot <- rbind(W.boot, W[first:last, ])
        }
      }
    } # for

    if (dims["MCD"] == 1) cov.method <- "mcd"
    if (dims["ACD"] == 1) cov.method <- "acd"
    if (dims["HPC"] == 1) cov.method <- "hpc"

    control <- jmcmControl()

    opt <-optimizeJmcm(m.boot, Y.boot, X.boot, Z.boot, W.boot,
      cov.method = cov.method, control = control, start = theta)

    result <- rbind(result, drop(opt$par))
    cat("iter ", iter, ": ", format(round(result[iter, ], 4), nsmall=4), "\n")

    beta.boot   <- drop(opt$par)[1:lbta]
    lambda.boot <- drop(opt$par)[(lbta + 1):(lbta + llmd)]
    gamma.boot  <- drop(opt$par)[(lbta + llmd + 1):(lbta + llmd + lgma)]

    Yest.boot <- rbind(Yest.boot, drop(X.ts %*% beta.boot))
    Zlmd.boot <- rbind(Zlmd.boot, drop(Z.ts %*% lambda.boot))
    Wgma.boot <- rbind(Wgma.boot, drop(W.tslag %*% gamma.boot))
  }

  Yest.boot <- apply(Yest.boot, 2, function(x) sort(x))
  Zlmd.boot <- apply(Zlmd.boot, 2, function(x) sort(x))
  Wgma.boot <- apply(Wgma.boot, 2, function(x) sort(x))

  Yest.u <- drop(Yest.boot[floor(0.975 * nboot), ])
  Yest.l <- drop(Yest.boot[ceiling(0.025 * nboot), ])

  Zlmd.u <- drop(Zlmd.boot[floor(0.975 * nboot), ])
  Zlmd.l <- drop(Zlmd.boot[ceiling(0.025 * nboot), ])
  Wgma.u <- drop(Wgma.boot[floor(0.975 * nboot), ])
  Wgma.l <- drop(Wgma.boot[ceiling(0.025 * nboot), ])

  plot(time, Y, xlab = "Time", ylab = "Response")
  lines(ts, Yest)
  lines(ts, Yest.u, lty = 2, lwd = 2)
  lines(ts, Yest.l, lty = 2, lwd = 2)

  if (dims["MCD"] == 1 || dims["ACD"] == 1) {
    xlab="Time"
    ylab="Log-innovat. var."
  }
  if (dims["HPC"] == 1) {
    xlab="Time"
    ylab="Log-variance"
  }
  plot(ts, Zlmd, type = 'l', xlab = xlab, ylab = ylab)
  lines(ts, Zlmd.u, lty = 2, lwd = 2)
  lines(ts, Zlmd.l, lty = 2, lwd = 2)

  if (dims["MCD"] == 1) {
    xlab="Lag"
    ylab="Autoregres. coeffic."
  }
  if (dims["ACD"] == 1) {
    xlab="Lag"
    ylab="MA. coeffic."
  }
  if (dims["HPC"] == 1) {
    xlab="Lag"
    ylab="Angles"
  }
  plot(tslag, Wgma, type = 'l', xlab = xlab, ylab = ylab)
  lines(tslag, Wgma.u, lty = 2, lwd = 2)
  lines(tslag, Wgma.l, lty = 2, lwd = 2)

}

# #' @export
# globalSearch <- function(formula, data = NULL,
#   cov.method = c('mcd', 'acd', 'hpc'),
#   control = jmcmControl())
# {
#   args <- ldFormula(formula, data, triple = c(1, 1, 1),
#                     cov.method, control)
#
#   m <- max(args$m)
#
#   triple <- c(m-1, m-1, m-1)
#   full <- ans <- jmcm(formula, data, triple=triple, cov.method, control)
#
#   bta0 <- ans@opt$beta
#   lmd0 <- ans@opt$lambda
#   gma0 <- ans@opt$gamma
#
#   cat("-------------------------------------------------------\n")
#
#   for(i in (m-2):1) {
#     triple <- c(i, m-1, m-1)
#     fit <- jmcm(formula, data, triple=triple, cov.method, control)
#     if(ans@opt$BIC > fit@opt$BIC) {
#       p   <- i
#       ans <- fit
#       cat("triple: ")    ; print(triple)
#       cat("    logLik: "); print(fit@opt$loglik)
#       cat("    BIC   : "); print(fit@opt$BIC)
#     }
#   }
#
#   cat("-------------------------------------------------------\n")
#
#   {
#     ans <- full
#     for(i in (m-2):1) {
#       triple <- c(m-1, i, m-1)
#       fit <- jmcm(formula, data, triple=triple, cov.method, control)
#       if(ans@opt$BIC > fit@opt$BIC) {
#         d <- i
#         ans <- fit
#         cat("triple: "); print(triple)
#         cat("    logLik: "); print(fit@opt$loglik)
#         cat("    BIC   : "); print(fit@opt$BIC)
#       }
#     }
#
#     cat("-------------------------------------------------------\n")
#
#     ans <- full
#     for(i in (m-2):1) {
#       triple <- c(m-1, m-1, i)
#       fit <- jmcm(formula, data, triple=triple, cov.method, control)
#       if(ans@opt$BIC > fit@opt$BIC) {
#         q <- i
#         ans <- fit
#         cat("triple: "); print(triple)
#         cat("    logLik: "); print(fit@opt$loglik)
#         cat("    BIC   : "); print(fit@opt$BIC)
#       }
#     }
#   }
#
#   cat("-------------------------------------------------------\n")
#
#   cat("p = "); print(p)
#   cat("d = "); print(d)
#   cat("q = "); print(q)
#
#   c(p,d,q)
# }
