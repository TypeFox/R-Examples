"geweke.diag" <-  function (x, frac1 = 0.1, frac2 = 0.5) 
{
    if (frac1 < 0 || frac1 > 1) {
        stop("frac1 invalid")
    }
    if (frac2 < 0 || frac2 > 1) {
        stop("frac2 invalid")
    }
    if (frac1 + frac2 > 1) {
        stop("start and end sequences are overlapping")
    }
    if (is.mcmc.list(x)) {
        return(lapply(x, geweke.diag, frac1, frac2))
    }
    x <- as.mcmc(x)
    xstart <- c(start(x), floor(end(x) - frac2 * (end(x) - start(x))))
    xend <- c(ceiling(start(x) + frac1 * (end(x) - start(x))), end(x))
    y.variance <- y.mean <- vector("list", 2)
    for (i in 1:2) {
        y <- window(x, start = xstart[i], end = xend[i])
        y.mean[[i]] <- apply(as.matrix(y), 2, mean)
        y.variance[[i]] <- spectrum0.ar(y)$spec/niter(y)
    }
    z <- (y.mean[[1]] - y.mean[[2]])/sqrt(y.variance[[1]] + y.variance[[2]])
    out <- list(z = z, frac = c(frac1, frac2))
    class(out) <- "geweke.diag"
    return(out)
}

"geweke.plot" <-
  function (x, frac1 = 0.1, frac2 = 0.5, nbins = 20, 
            pvalue = 0.05, auto.layout = TRUE, ask, ...) 
{
  if (missing(ask)) {
    ask <- if (is.R()) {
      dev.interactive()
    }
    else {
      interactive()
    }
  }
    
  x <- as.mcmc.list(x)
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) 
    oldpar <- par(mfrow = set.mfrow(Nchains = nchain(x), 
                    Nparms = nvar(x)))
  ystart <- seq(from = start(x), to = (start(x) + end(x))/2, length = nbins)
  if (is.R())
    gcd <- array(dim = c(length(ystart), nvar(x), nchain(x)), 
               dimnames = c(ystart, varnames(x), chanames(x)))
  else
    gcd <- array(dim = c(length(ystart), nvar(x), nchain(x)), 
               dimnames = list(ystart, varnames(x), chanames(x)))
               
  for (n in 1:length(ystart)) {
    geweke.out <- geweke.diag(window(x, start = ystart[n]), 
                              frac1 = frac1, frac2 = frac2)
    for (k in 1:nchain(x)) gcd[n, , k] <- geweke.out[[k]]$z
  }
  climit <- qnorm(1 - pvalue/2)
  for (k in 1:nchain(x)) for (j in 1:nvar(x)) {
    ylimit <- max(c(climit, abs(gcd[, j, k])))
    plot(ystart, gcd[, j, k], type = "p", xlab = "First iteration in segment", 
         ylab = "Z-score", pch = 4, ylim = c(-ylimit, ylimit), 
         ...)
    abline(h = c(climit, -climit), lty = 2)
    if (nchain(x) > 1) {
      title(main = paste(varnames(x, allow.null = FALSE)[j], 
              " (", chanames(x, allow.null = FALSE)[k], ")", 
              sep = ""))
    }
    else {
      title(main = paste(varnames(x, allow.null = FALSE)[j], 
              sep = ""))
    }
    if (k==1 && j==1)
       oldpar <- c(oldpar, par(ask = ask))
  }
  invisible(list(start.iter = ystart, z = gcd))
}

"print.geweke.diag" <- function (x, digits = min(4, .Options$digits), ...) 
  ## Print method for output from geweke.diag
{
  cat("\nFraction in 1st window =", x$frac[1])
  cat("\nFraction in 2nd window =", x$frac[2], "\n\n")
  print.default(x$z, digits = digits, ...)
  cat("\n")
  invisible(x)
}











