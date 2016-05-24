svsim <- function(len, mu = -10, phi = 0.98, sigma = 0.2, nu = Inf) {
 
 # Some error checking
 if (any(is.na(len)) | !is.numeric(len) | length(len) != 1 | any(len < 1)) {
  stop("Argument 'len' (length of simulated series) must be a single number >= 2.")
 } else {
  len <- as.integer(len)
 }

 if (!is.numeric(mu) | length(mu) != 1) {
  stop("Argument 'mu' (level of latent variable) must be a single number.")
 }

if (!is.numeric(phi) | length(phi) != 1) {
  stop("Argument 'phi' (persistence of latent variable) must be a single number.")
 }

if (!is.numeric(sigma) | length(sigma) != 1 | sigma <= 0) {
  stop("Argument 'sigma' (volatility of latent variable) must be a single number > 0.")
 }

if (!is.numeric(nu) || length(nu) != 1 || nu <= 2) {
 stop("Argument 'nu' (degrees of freedom for the conditional error) must be a single number > 2.")
}

 h <- rep(as.numeric(NA), len)
 h0 <- rnorm(1, mean=mu, sd=sigma/sqrt(1-phi^2))
 innov <- rnorm(len)

 # simulate w/ simple loop
 h[1] <- mu + phi*(h0-mu) + sigma*innov[1]
 for (i in seq(2, len = len-1)) h[i] <- mu + phi*(h[i-1]-mu) + sigma*innov[i]

 if (is.infinite(nu)) {
  y <- exp(h / 2) * rnorm(len)  # "log-returns"
 } else {
  y <- exp(h / 2) * rt(len, df = nu)  # "log-returns"
 }
 ret <- list(y = y, vol = exp(h/2), vol0 = exp(h0/2),
	     para = list(mu = mu,
			 phi = phi,
			 sigma = sigma))
 if (is.finite(nu)) ret$para$nu <- nu
 class(ret) <- "svsim"
 ret
}

print.svsim <- function(x, ...) {
 cat("\nSimulated time series consisting of", length(x$y), "observations.\n
Parameters: level of latent variable                  mu =", x$para$mu, "
            persistence of latent variable           phi =", x$para$phi, "
            standard deviation of latent variable  sigma =", x$para$sigma, "
            ")
 if (length(x$para) == 4) cat("degrees of freedom parameter              nu =", x$para$nu, "
	    ")
 cat("\nSimulated initial conditional volatility:", x$vol0, "\n")
 cat("\nSimulated conditional volatilities:\n")
 print(x$vol, ...)
 cat("\nSimulated data (usually interpreted as 'log-returns'):\n")
 print(x$y, ...)
}

plot.svsim <- function(x, mar = c(3, 2, 2, 1), mgp = c(1.8, .6, 0), ...) {
 op <- par(mfrow = c(2, 1), mar = mar, mgp = mgp)
 plot.ts(100*x$y, ylab = "", ...)
 mtext("Simulated data: 'log-returns' (in %)", cex = 1.2, line = .4, font = 2)
 plot.ts(100*x$vol, ylab = "", ...)
 mtext("Simulated conditional volatilities (in %)", cex = 1.2, line = .4, font = 2)
 par(op)
}

summary.svsim <- function(object, ...) {
 ret <- vector("list")
 class(ret) <- "summary.svsim"
 ret$len <- length(object$y)
 ret$para <- object$para
 ret$vol0 <- 100*object$vol0
 ret$vol <- summary(100*object$vol)
 ret$y <- summary(100*object$y)
 ret
}

print.summary.svsim  <- function(x, ...) {
 cat("\nSimulated time series consisting of ", x$len, " observations.\n",
     "\nParameters: level of latent variable                  mu = ",
     x$para$mu, 
     "\n            persistence of latent variable           phi = ",
     x$para$phi,
     "\n            standard deviation of latent variable  sigma = ",
     x$para$sigma, "\n", sep="")
 if (length(x$para) == 4) cat("            degrees of freedom parameter              nu =", x$para$nu, "
	    ")

 cat("\nSimulated initial conditional volatility (in %): ")
 cat(x$vol0, "\n")
 cat("\nSummary of simulated conditional volatilities (in %):\n")
 print(x$vol)
 cat("\nSummary of simulated data (in %):\n")
 print(x$y)
 invisible(x)
}
