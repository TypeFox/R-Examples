JAGSrun <-  function(y, prefix = yname, model = BMMmodel(k = 2),
                     control = JAGScontrol(variables = c("mu", "tau", "eta")),
                     tmp = TRUE, cleanup = TRUE, ...) {
  yname <- deparse(substitute(y))
  if (!is.null(dim(y))) {
    if (dim(y)[1] == 1) y <- y[1,]
    else if (dim(y)[2] == 1) y <- y[,1]
    else stop("Only univariate response allowed")
  }
  y <- as.numeric(y)
  cl <- match.call()
  if (tmp) {
    dir <- getwd()
    tmpdir <- tempdir()
    if (!file.exists(tmpdir)){
      if (!dir.create(tmpdir)) stop("Error creating tmp directory")
    }
    setwd(tmpdir)
  }
  results <- JAGScall(model, y, prefix, control, ...)
  if (cleanup) {
    unlink(paste(prefix, ".bug", sep = ""))
  }
  if (tmp) setwd(dir)
  z <- list(call = cl, results = results$results, model = results$model,
            variables = results$variables, data = y)
  class(z) <- "jags"
  z
}

# Print method for jags objects adapted from print.mcmc and
# summary.mcmc in package coda written by Martyn Plummer, Nicky Best,
# Kate Cowles, Karen Vines

print.jags <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (inherits(x$model, "BMMmodel")) {
    if (is.null(x$results)) cat("No results!\n")
    else {
      cat("Markov Chain Monte Carlo (MCMC) output:\nStart =", stats::start(x$results), 
          "\nEnd =", stats::end(x$results), "\nThinning interval =", thin(x$results), 
          "\n")
      for (i in x$variables) {
        y <- x$results[,grep(paste("^", i, sep = ""), colnames(x$results)), drop = FALSE]
        if(dim(y)[2] <=  x$model$data$k) {
          yout <- summaryShort.mcmc(y)
          class(yout) <- "summaryShort.mcmc"
          cat(paste("\n Empirical mean, standard deviation and 95% CI for", i, "\n"))
          print(yout, ...)
        }
      }
    }
  }
}

summaryShort.mcmc <- function (object, quantiles = c(0.025, 0.975), 
                               ...) {
  x <- methods::as(object, "mcmc")
  statnames <- c("Mean", "SD")
  varstats <- matrix(nrow = nvar(x), ncol = length(statnames), 
                     dimnames = list(varnames(x), statnames))
  if (is.matrix(x)) {
    xmean <- apply(x, 2, mean)
    xvar <- apply(x, 2, stats::var)
    varquant <- t(apply(x, 2, stats::quantile, quantiles))
  }
  else {
    xmean <- mean(x, na.rm = TRUE)
    xvar <- stats::var(x, na.rm = TRUE)
    varquant <- stats::quantile(x, quantiles)
  }
  varstats[, 1] <- xmean
  varstats[, 2] <- sqrt(xvar)
  varstats <- drop(varstats)
  varquant <- drop(varquant)
  out <- list(statistics = varstats, quantiles = varquant,
              start = stats::start(x), end = stats::end(x), thin = thin(x), nchain = 1)
  class(out) <- "summaryShort.mcmc"
  return(out)
}

print.summaryShort.mcmc <- function(x, digits = max(3, .Options$digits - 3), ...) {
  if (is.matrix(x$statistics)) {
    print(cbind(x$statistics, x$quantiles), digits = digits, ...)
  }
  else print(c(x$statistics, x$quantiles), digits = digits, ...)
}

