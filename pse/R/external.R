#########################################
############# External code #############
#########################################

# The functions in this file were copied from other
# packages, either to reduce the number of package
# dependencies, or because they are not exported in their
# package namespace. Original authors and licenses are
# given for all of them.

# Nodeplot by Gilles Pujol 2006
# Released in package sensitivity as GPL-2
nodeplot <-
function(x, xlim = NULL, ylim = NULL, labels = TRUE,
                     col = par("col"), pch = 22, bg = "white",
                     add = FALSE, at = NULL, ...) {
  n <- nrow(x)
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (add) {
    par(new = TRUE)
  }

  # axes
  
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", type = "n", ...)
  if (class(labels) == "logical") {
    if (labels) {
      axis(side = 1, at = at, labels = rownames(x))
    } else {
      axis(side = 1, at = at, labels = FALSE, tick = FALSE)
    }
  } else if (class(labels) == "character") {
    axis(side = 1, at = at, labels = labels)
  }
  axis(side = 2)
  box()

  # bias

  if ("bias" %in% colnames(x)) {
	  # added by A.C. 2014: if bias is NA, it should not be added
	  x[["bias"]][is.na(x[["bias"]])]  <- 0
    xx <- x[["original"]] - x[["bias"]]
  } else {
    xx <- x[["original"]]
  }
  
  # confidence intervals

  if (("min. c.i." %in% colnames(x)) & "max. c.i." %in% colnames(x)) {
    for (i in 1 : n) {
      lines(c(at[i], at[i]), c(x[["min. c.i."]][i], x[["max. c.i."]][i]),
            col = col)
    }
  }

  # points

  points(at, xx, col = col, pch = pch, bg = bg)
}

# Vec operator by Tarn Duong
# Released in package ks, licensed as GPL-2/3
vec <- function(x, byrow=FALSE)
{
	  if (is.vector(x)) return (x)

  if (byrow) x <- t(x)
    d <- ncol(x)
    vecx <- vector()
	  for (j in 1:d)
		      vecx <- c(vecx, x[,j])

	  return(vecx)
}

#        Bootstrap statistics (overlay for the boot package)
#                         Gilles Pujol 2006
#        Released in package sensitivity, licensed as GPL-2
# Edited to make the NAMESPACE less bloated

# bootstats(b, conf = 0.95, type = "norm")
# b : object of class 'boot'
# confidence : confidence level for bootstrap bias-corrected confidence
#   intervals
# type : type of confidence interval, "norm" or "basic"
#
# returns : a data.frame of bootstrap statistics

bootstats <- function(b, conf = 0.95, type = "norm") {
  p <- length(b$t0)
  lab <- c("original", "bias", "std. error", "min. c.i.", "max. c.i.")
  out <-  as.data.frame(matrix(nrow = p, ncol = length(lab),
                               dimnames = list(NULL, lab)))

  for (i in 1 : p) {
    
    # original estimation, bias, standard deviation
      
    out[i, "original"] <- b$t0[i]
    out[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    out[i, "std. error"] <- sd(b$t[, i])
      
    # confidence interval

    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)
      if (!is.null(ci)) {
        out[i, "min. c.i."] <- ci$norm[2]
        out[i, "max. c.i."] <- ci$norm[3]
      }
    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        out[i, "min. c.i."] <- ci$basic[4]
        out[i, "max. c.i."] <- ci$basic[5]
      }
    }
  }
  
  return(out)
}
# Partial Correlation Coefficients
# Released originally as package sensitivity, licensed as GPL-2
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

# Based on code by Gilles Pujol, modified by Andre Chalom 2013
pcc <- function(X, y, rank, nboot, conf, ...) UseMethod("pcc")
pcc.LHS <- function(X, y=NULL, rank=FALSE, nboot=0, conf=0.95, ...) {
	res <- get.results(X)
	L <- get.data(X)
	f <- function(r) pcc(X=L, y=r, rank, nboot, conf, ...)
	apply(res, 2, f)
}
pcc.default <- function(X, y=NULL, rank = FALSE, nboot = 0, conf = 0.95, ...) {
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
    boot.pcc <- boot::boot(data, estim.pcc, R = nboot)
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

#' External functions
#' 
#' This function is derived from package "sensitivity", please proceed to its help files.
#' @param x A PCC or PRCC object
#' @param \dots Currently ignored
#' @export
#' @rdname external
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

# Based on code by Gilles Pujol 2006, modified by Andre Chalom 2013
plot.pcc <-
function(x, ylim = c(-1,1), main=NULL, ...) {  
		main = (if ("PCC" %in% names(x)) "PCC" else "PRCC")
		if ("PCC" %in% names(x)) {
				nodeplot(x$PCC, ylim = ylim, main=main, ...)
		}else if ("PRCC" %in% names(x)) {
				nodeplot(x$PRCC, ylim = ylim, main=main, ...)
		}  
}
