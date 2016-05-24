# Standardized Regression Coefficients
#
# Gilles Pujol 2006


estim.src <- function(data, i = 1:nrow(data)) {
  d <- data[i, ]
  lm.Y <- lm(formula(paste(colnames(d)[1], "~", paste(colnames(d)[-1], collapse = "+"))), data = d)
  coefficients(lm.Y)[-1] * sapply(d[-1], sd) / sapply(d[1], sd)
}


src <- function(X, y, rank = FALSE, nboot = 0, conf = 0.95) {
  data <- data.frame(Y = y, X)

  if (rank) {
    for (i in 1:ncol(data)) {
      data[,i] <- rank(data[,i])
    }
  }
  
  if (nboot == 0) {
    src <- data.frame(original = estim.src(data))
    rownames(src) <- colnames(X)
  } else {
    boot.src <- boot(data, estim.src, R = nboot)
    src <- bootstats(boot.src, conf, "basic")
    rownames(src) <- colnames(X)
  }

  out <- list(X = X, y = y, rank = rank, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "src"
  if (! rank) {
    out$SRC <- src
  } else {
    out$SRRC = src
  }
  return(out)
}


print.src <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if ("SRC" %in% names(x)) {
    cat("\nStandardized Regression Coefficients (SRC):\n")
    print(x$SRC)
  } else if ("SRRC" %in% names(x)) {
    cat("\nStandardized Rank Regression Coefficients (SRRC):\n")
    print(x$SRRC)
  }
}


plot.src <- function(x, ylim = c(-1,1), ...) {  
  if ("SRC" %in% names(x)) {
    nodeplot(x$SRC, ylim = ylim, main = "SRC")
  } else if ("SRRC" %in% names(x)) {
    nodeplot(x$SRRC, ylim = ylim, main = "SRRC")
  }
}
