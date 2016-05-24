# Partial Correlation Coefficients
# Gilles Pujol 2006


estim.pcc <- function(data, i = 1:nrow(data)) {  
  d <- data[i, ]
  p <- ncol(d) - 1
  pcc <- numeric(p)
  for (j in 1:p) {
    Xtildej.lab <- paste(colnames(d)[c(-1, -j-1)], collapse = "+")
    lm.Y <- lm(formula(paste(colnames(d)[1], "~", Xtildej.lab)), data = d)
    lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
    pcc[j] <- cor(d[1] - fitted(lm.Y), d[j+1] - fitted(lm.Xj))
  }
  pcc
}


pcc <- function(X, y, rank = FALSE, nboot = 0, conf = 0.95) {
  data <- cbind(Y = y, X)

  if (rank) {
    for (i in 1:ncol(data)) {
      data[,i] <- rank(data[,i])
    }
  }
  
  if (nboot == 0) {
    pcc <- data.frame(original = estim.pcc(data))
    rownames(pcc) <- colnames(X)
  } else {
    boot.pcc <- boot(data, estim.pcc, R = nboot)
    pcc <- bootstats(boot.pcc, conf, "basic")
    rownames(pcc) <- colnames(X)
  }

  out <- list(X = X, y = y, rank = rank, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "pcc"
  if (! rank) {
    out$PCC <- pcc
  } else {
    out$PRCC = pcc
  }
  return(out)
}


print.pcc <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if ("PCC" %in% names(x)) {
    cat("\nPartial Correlation Coefficients (PCC):\n")
    print(x$PCC)
  } else if ("PRCC" %in% names(x)) {
    cat("\nPartial Rank Correlation Coefficients (PRCC):\n")
    print(x$PRCC)
  }
}


plot.pcc <- function(x, ylim = c(-1,1), ...) {  
  if ("PCC" %in% names(x)) {
    nodeplot(x$PCC, ylim = ylim, main = "PCC")
  }else if ("PRCC" %in% names(x)) {
    nodeplot(x$PRCC, ylim = ylim, main = "PRRC")
  }
}
