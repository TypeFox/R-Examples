# Sobol' indices estimation (Sobol 1993)
#
# Gilles Pujol 2006


# subsets : subsets of size (less or) equal to n

subsets <- function(set, size = length(set), type = "leq") {
  if (size <= 0) {
    NULL
  } else if (size == 1) {
    as.list(set)
  } else if (type == "leq") {
    c(Recall(set, size - 1, "leq"), Recall(set, size, "eq"))
  } else if (size == length(set)) {
    list(set)
  } else {
    c(lapply(Recall(set[-1], size - 1, "eq"), function(x) c(set[1], x)),
      Recall(set[-1], size, "eq"))
  }
}


sobol <- function(model = NULL, X1, X2, order = 1, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) {
    stop("The samples X1 and X2 must have the same dimensions")
  }
  p <- ncol(X1)

  if (is.numeric(order) && (length(order) == 1)) {
    indices.list <- subsets(set = 1 : p, size = order)
  } else if (is.list(order)) {
    indices.list <- order
  } else {
    stop("The argument \'order\' must be either a scalar or a list")
  }
  
  X <- X1
  for (i in indices.list) {
    Xb <- X2
    Xb[, i] <- X1[, i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, order = order, nboot = nboot,
            conf = conf, X = X, call = match.call())
  class(x) <- "sobol"
  
  if (! is.null(x$model)) {
    response(x, ...)
    tell(x, ...)
  }
  
  return(x)
}


estim.sobol <- function(data, i = 1 : nrow(data)) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  V <- var(d[, 1])
  VCE <- colSums(d[, -1] * d[, 1]) / (n - 1) - mean(d[, 1])^2
  c(V, VCE)
}


tell.sobol <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))

  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }

  p <- ncol(x$X1)
  n <- nrow(x$X1)
  
  if (is.numeric(x$order)) {
    indices.list <- subsets(set = 1 : p, size = x$order)
  } else if (is.list(x$order)) {
    indices.list <- x$order
  }
  ni <- length(indices.list)
  
  indices.labels <- lapply(indices.list, function(i) paste(colnames(x$X1)[i], collapse="*"))
  
  data <- matrix(x$y, nrow = n)

  # estimation of the variances of the conditional expectations (V)

  if (x$nboot == 0) {
    V <- data.frame(original = estim.sobol(data))
  } else {
    V.boot <- boot(data, estim.sobol, R = x$nboot)
    V <- bootstats(V.boot, x$conf, "basic")
  }
  rownames(V) <- c("global", indices.labels)

  # estimation of the terms of the ANOVA decomposition (D)
  # and the Sobol' indices (S)
  
  if (is.numeric(x$order)) {

    # 'mask' is a matrix containing 0, +1 and -1 and is used to compute the
    # terms of the ANOVA decomposition. (Indeed, D_I = V_I - sum D_J for J
    # subset of I, and #J < # I ; +1 correspond to the set I and -1 correspond
    # to the subsets Js). 
    mask <- matrix(0, nrow = ni, ncol = ni)
    for (i in 1 : ni) {
      mask[i, i] <- 1
      lower.i <- subsets(set = indices.list[[i]],
                         size = length(indices.list[[i]]) - 1)
      mask[i, match(lower.i, indices.list)] <- -1
    }

    if (x$nboot == 0) {
      D <- V[-1, 1, drop = FALSE]
      for (i in 1:ni) {
        D[i,1] <- sum(D * mask[i,])
      }

      S <- D / V[1,1]
    } else {
      D.boot <- V.boot
      D.boot$t0 <- V.boot$t0[-1]
      D.boot$t <- V.boot$t[,-1, drop = FALSE]
      for (i in 1:ni) {
        D.boot$t0[i] <- sum(D.boot$t0 * mask[i,])
        D.boot$t[,i] <- colSums(t(D.boot$t) * mask[i,])
      }
      D <- bootstats(D.boot, x$conf, "basic")
      rownames(D) <- indices.labels

      S.boot <- D.boot
      S.boot$t0 <- D.boot$t0 / V.boot$t0[1]
      S.boot$t <- D.boot$t / V.boot$t[,1]
      S <- bootstats(S.boot, x$conf, "basic")
      rownames(S) <- indices.labels
    }
  }

  # return
  x$V <- V
  if (is.numeric(x$order)) {
    x$D <- D
    x$S <- S
  }

  for (i in return.var) {
    x[[i]] <- get(i)
  }

  assign(id, x, parent.frame())
}


print.sobol <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      cat("\nSobol indices\n")
      print(x$S)
    }
  } else {
    cat("(empty)\n")
  }
}


plot.sobol <- function(x, ylim = c(0, 1), ...) {
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
  }
}
