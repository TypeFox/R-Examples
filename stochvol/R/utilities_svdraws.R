logret <- function(x, demean = FALSE) {
 tmp <- diff(log(x))
 if (is(tmp, "xts")) tmp <- tmp[-1]
 if (all(isTRUE(demean))) tmp <- tmp - mean(tmp)
 tmp
}

para <- function(x) {
 x$para
}

latent <- function(x) {
 x$latent
}

latent0 <- function(x) {
 x$latent0
}

priors <- function(x) {
 x$priors
}

thinning <- function(x) {
 x$thinning
}

runtime <- function(x) {
 x$runtime
}

updatesummary <- function(x, quantiles = c(.05, .5, .95), esspara = TRUE, esslatent = FALSE) {
 
 # NEW: Check if conditional t errors are used
 if (ncol(x$para) == 4) terr <- TRUE else terr <- FALSE
 if (exists("betas", x)) regression <- TRUE else regression <- FALSE

 summaryfunction <- function(x, quants = quantiles, ess = TRUE) {
  if (ess) {
   c(mean = mean(x), sd = sd(x), quantile(x, quantiles),
    ESS = as.numeric(effectiveSize(x)))
  } else {
   c(mean = mean(x), sd = sd(x), quantile(x, quantiles))
  }
 }

 res <- list()
 
 res$para <- t(apply(x$para, 2, summaryfunction, ess = esspara))
 res$para <- rbind(res$para, "exp(mu/2)" = c(summaryfunction(exp(x$para[,"mu"]/2), ess=FALSE), res$para["mu", "ESS"]))
 res$para <- rbind(res$para, "sigma^2" = c(summaryfunction(x$para[,"sigma"]^2, ess=FALSE), res$para["sigma", "ESS"]))
 
 res$latent <- t(apply(x$latent, 2, summaryfunction, ess = esslatent))
 tmp <- exp(x$latent/2)
 res$latent <- cbind(res$latent, "mean(exp(h_t/2))" = colMeans(tmp))
 res$latent <- cbind(res$latent, "sd(exp(h_t/2))" = apply(tmp, 2, sd))
 
 res$latent0 <- c(summaryfunction(x$latent0, ess = esslatent), "mean(exp(h_t/2))" = mean(exp(x$latent0/2)), "sd(exp(h_t/2))" = sd(exp(x$latent0/2)))

 if (terr && x$thinning$para == x$thinning$latent) {
  tmp <- sqrt((exp(x$latent) * x$para[,"nu"]) / (x$para[,"nu"] - 2))
  res$sd <- t(apply(tmp, 2, summaryfunction, ess = esslatent))
  rownames(res$sd) <- gsub("h", "sd", rownames(res$latent))
 }

 if (exists("beta", x)) res$beta <- t(apply(x$beta, 2, summaryfunction, ess = esspara))
 
 x$summary <- res
 x
}

summary.svdraws <- function(object, showpara = TRUE, showlatent = TRUE, ...) {
 ret <- vector("list")
 class(ret) <- "summary.svdraws"
 ret$mcp <- mcpar(para(object))
 ret$mcl <- mcpar(latent(object))
 ret$priors <- priors(object)
 if (all(isTRUE(showpara))) ret$para <- para(object$summary)
 if (all(isTRUE(showlatent))) ret$latent <- rbind("h_0" = latent0(object$summary), latent(object$summary))
 ret
}

print.summary.svdraws <- function(x, ...) { 
 cat("\nSummary of ", x$mcp[2]-x$mcp[1]+x$mcp[3], " MCMC draws after a burn-in of ", x$mcp[1]-x$mcp[3], ".\n
Prior distributions:
mu        ~ Normal(mean = ", x$priors$mu[1], ", sd = ", x$priors$mu[2], ")
(phi+1)/2 ~ Beta(a0 = ", x$priors$phi[1], ", b0 = ", x$priors$phi[2], ")
sigma^2   ~ ", x$priors$sigma, " * Chisq(df = 1)\n", sep="")

 if (length(x$priors) == 4) cat("nu        ~ Unif(lower = ", x$priors$nu[1], ", upper = ", x$priors$nu[2], ")\n", sep="")

 if (exists("para", x)) {
  cat("\nPosterior draws of parameters (thinning = ", x$mcp[3], "):\n", sep='')
  print(x$para, digits=2, ...)
 }
 
 if (exists("latent", x)) {
  cat("\nPosterior draws of initial and contemporaneous latents (thinning = ", x$mcl[3], "):\n", sep='')
  print(x$latent, digits=2, ...)
 }
 invisible(x)
}

print.svdraws <- function(x, showpara = TRUE, showlatent = TRUE, ...) {
 if (all(isTRUE(showpara))) {
  cat("\n*** Posterior draws of parameters ***\n")
  print(para(x), ...)
 }
 
 if (all(isTRUE(showlatent))) {
  cat("\n*** Posterior draws of initial latent variable h_0 ***\n")
  print(latent0(x), ...)
  cat("\n*** Posterior draws of contemporaneous latent variables h_t ***\n")
  print(latent(x), ...)
 }
 invisible(x)
}

# residuals
residuals.svdraws <- function(object, type = "mean", ...) {
 if (!is(object, "svdraws")) stop("This function expects an 'svdraws' object.")
 if (!any(type == "mean", type == "median")) stop("Argument 'type' must currently be either 'mean' or 'median'.")

 if (object$thinning$time != 1) warning("Not every point in time has been stored ('thintime' was set to a value unequal to 1 during sampling), thus only some residuals have been extracted.")
  
 if (type == "mean") {
  res <- rowMeans(as.numeric(object$y)[seq(1, length(object$y), by=object$thinning$time)] / exp(t(object$latent)/2))
 }
 
 if (type == "median") {
  res <- apply(as.numeric(object$y)[seq(1, length(object$y), by=object$thinning$time)] / exp(t(object$latent)/2), 1, median)
 }

 names(res) <- sub("h", "r", colnames(object$latent))
 class(res) <- "svresid"
 attr(res, "type") <- type

 # NEW: Also return posterior mean/median of df parameter if terr = TRUE
 if (ncol(object$para) == 4) attr(res, "nu") <- get(type)(object$para[,"nu"])

 res
}

predict.svdraws <- function(object, steps = 1L, ...) {
 if (!is(object, "svdraws")) stop("Argument 'object' must be of class 'svdraws'.")
 steps <- as.integer(steps)
 if (steps < 1) stop("Argument 'steps' must be greater or equal to 1.")
 thinlatent <- object$thinning$latent
 thinpara <- object$thinning$para
 if (thinpara != thinlatent) {
  warning("Thinning of parameters is different from thinning of latent variables. Trying to sort this out.")
  if (thinpara %% thinlatent == 0) {
   usepara <- 1:(dim(object$para)[1])
   uselatent <- seq(thinpara/thinlatent, dim(object$latent)[1], by=thinpara/thinlatent)
  } else if (thinlatent %% thinpara == 0) {
   uselatent <- 1:(dim(object$latent)[1])
   usepara <- seq(thinlatent/thinpara, dim(object$para)[1], by=thinlatent/thinpara)
  } else stop("Incompatible thinning parameters. Prediction currently not implemented.")
 } else {
  usepara <- uselatent <- seq.int(dim(object$para)[1])
 }
 
 mu <- object$para[,"mu"][usepara]
 phi <- object$para[,"phi"][usepara]
 sigma <- object$para[,"sigma"][usepara]
 hlast <- object$latent[,dim(object$latent)[2]][uselatent]

 mythin <- max(thinpara, thinlatent)
 len <- length(sigma)
 volpred <- mcmc(matrix(as.numeric(NA), nrow=len, ncol=steps), start=mythin, end=len*mythin, thin=mythin)
 
 volpred[,1] <- mu + phi*(hlast - mu) + sigma*rnorm(len)
 if (steps > 1)
  for (i in (seq.int(steps-1) + 1))
   volpred[,i] <- mu + phi*(volpred[,i-1] - mu) + sigma*rnorm(len)

 class(volpred) <- c("svpredict", "mcmc")
 lastname <- dimnames(object$latent)[[2]][dim(object$latent)[2]]
 lastnumber <- as.integer(gsub("h_", "", lastname))
 colnames(volpred) <- paste("h_", seq(lastnumber + 1, lastnumber + steps), sep='')
 volpred
}

# used to forecast AR-SV models (needs more testing!)
arpredict <- function(object, volpred) {
 if (!is(object, "svdraws")) stop("Argument 'object' must be of class 'svdraws'.")
 if (!is(volpred, "svpredict")) stop("Argument 'volpred' must be of class 'svpredict'.")
 if (colnames(object$priors$designmatrix)[1] == "const") dynamic <- TRUE else stop("Probably not an AR-specification.")
 order <- ncol(object$priors$designmatrix) - 1

 len <- nrow(volpred)
 steps <- ncol(volpred)
 fromtoby <- attr(volpred, "mcpar")
 usepara <- seq(from = fromtoby[1], to = fromtoby[2], by = fromtoby[3])

 if (ncol(object$para) == 4) {
  nu <- object$para[,"nu"][usepara]
 } else {
  nu <- Inf  # corresponds to conditional normality
 }

 if (ncol(object$beta) > 1) {
  betarev <- object$beta[usepara,c(1,ncol(object$beta):2)]
  lastX <- matrix(c(1, object$y[length(object$y) - order:1 + 1]), nrow = 1)
 } else {
  betarev <- object$beta[usepara,,drop=FALSE]
  lastX <- matrix(1, nrow = 1)
 }

 meanpred <- mcmc(matrix(as.numeric(NA), nrow = len, ncol = steps),
		  start = fromtoby[1], end = fromtoby[2], thin = fromtoby[3])
 meanpred[,1] <- tcrossprod(lastX, betarev) + exp(volpred[,1]/2)*rt(len, df = nu)
 if (steps > 1) {
  lastX <- matrix(rep(lastX, len), nrow = len, byrow = TRUE)
  for (i in (seq.int(steps-1) + 1)) {
   if (ncol(object$beta) > 1)
    lastX <- cbind(lastX[,-2, drop = FALSE], meanpred[,i-1])
   meanpred[,i] <- rowSums(lastX*betarev) + exp(volpred[,i]/2)*rt(len, df = nu)
  }
 }
 class(meanpred) <- c("distpredict", "mcmc")
 meanpred
}
