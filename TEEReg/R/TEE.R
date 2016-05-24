TEE <- function (formula, data, offset = NULL, p.trimmed = NULL, p.subsample = 1, method = "tee") {
  # Error checks
  if (missing(formula)) {
    stop("'formula' must be provided.")
  }
  if (missing(data)) {
    stop("'data' must be provided.")
  }
  if (method != "ols" & method != "tee") {
     stop(gettextf("invalid 'method' argument, method = '%s' is not supported. Using 'tee' or 'ols'.", method), domain = NA)
  }
  if (is.null(p.trimmed) & method == "tee") {
     stop("'p.trimmed' must be provided when 'method' is 'tee'.")
  }
  if (!is.null(p.trimmed)){
     if (!is.numeric(p.trimmed)) {
       stop("'p.trimmed' must be numeric.")
     } else if (p.trimmed >= 1 | p.trimmed < 0) {
       stop("invalid 'p.trimmed' argument.")
     }
  }
  if (!is.numeric(p.subsample)) {
     stop("'p.subsample' must be numeric.")
  } else if (p.subsample > 1 | p.subsample <= 0) {
     stop("invalid 'p.subsample' argument.")
  }
  mcall <- match.call(expand.dots = FALSE) # returns a call in which all of the specified arguments are specified by their full names
  mat <- match(c("formula", "data", "offset"), names(mcall), 0L)  # returns a vector of the positions matches of its first argument in its second
  mcall <- mcall[c(1L, mat)]
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- quote(stats::model.frame)
  mcall <- eval(mcall, parent.frame())  # evaluate an R expression in a specified enviroment
  mcallt <- attr(mcall, "terms")
  if (!is.null(offset)) {
    if (length(offset) != nrow(data)) {
      stop(gettextf("number of offsets is %d, should equal %d (number of observations).", length(offset), nrow(data)), domain = NA)
    } else {
      offset <- as.vector(model.offset(mcall))
    }
  }
  Yall <- model.response(mcall, "any")
  # model does not contain intercept and covariates 
  if (is.empty.model(mcallt)) {
    Xall <- NULL
    output <- list(coefficients = if (is.matrix(Yall)) matrix(, 0, 3) else numeric(), residuals = Yall, fitted.values = 0*Yall, rank = 0L)
    if (is.null(offset)) {
      output$fitted.values <- offset
      output$residuals <- Yall - offset
    }
    print(list(output))
    stop("no parameters need to be estimated.") 
  } else {
    Xall <- model.matrix(mcallt, mcall)
    names <- colnames(Xall)
  }
  if (method == "ols") {
    callt <- match.call()
    c <- match(c("formula", "data", "offset", "method"), names(callt), 0L)  # returns a vector of the positions matches of its first argument in its second
    callt <- callt[c(1L, c)]
    if (abs(det(t(Xall)%*%Xall)) < 1e-08) {
      warning("Matrix is singular, generalized inverse is used.")
      tol = sqrt(.Machine$double.eps)
      XpXsvd <- svd(t(Xall)%*%Xall)
      Positive <- XpXsvd$d > max(tol * XpXsvd$d[1L], 0)
      if (all(Positive)) { 
        hat <- XpXsvd$v %*% (1/XpXsvd$d * t(XpXsvd$u))
      } else if (!any(Positive)) { 
        hat <- array(0, dim(Xall)[2L:1L])
      } else {
        hat <- XpXsvd$v[, Positive, drop = FALSE] %*% ((1/XpXsvd$d[Positive])*t(XpXsvd$u[, Positive, drop = FALSE]))
      }
    } else {
      hat <- solve(t(Xall)%*%Xall)
    }
    if (!is.null(offset)) {
      TEE.est <- as.matrix(t(hat%*%t(Xall)%*%(Yall-offset)))
    } else { 
      TEE.est <- as.matrix(t(hat%*%t(Xall)%*%Yall))
    }
  } else if (method == "tee") {
    samplesize <- length(Yall)
    p <- ncol(Xall)
    index <- combn((1:samplesize), p)
    k <- ncol(index)
    set.seed(23211342)
    s <- ceiling(p.subsample*k)
    subset <- as.matrix(index[,sample(1:k, s, replace = FALSE)])
    beta.h <- matrix(NA, nrow = p, ncol = s)
    det.XhXh <- c()
    sum.abse <- c()
    r <- round((1-p.trimmed)*s)
    for (i in 1:s) {
      Y <- Yall[subset[,i]]
      X <- Xall[subset[,i],]
      # Solve sigular problem using generalized inverse
      if (abs(det(X)) < 1e-08) {
        tol = sqrt(.Machine$double.eps)
        Xsvd <- svd(X)
        Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
        if (all(Positive)) { 
          hat <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
        } else if (!any(Positive)) { 
          hat <- array(0, dim(X)[2L:1L])
        } else {
          hat <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive])*t(Xsvd$u[, Positive, drop = FALSE]))
        }
      } else {
        hat <- solve(X)
      }
      if (!is.null(offset)) {
        beta.h[,i] <- hat %*% (Y-offset[subset[,i]])
      } else { 
        beta.h[,i] <- hat%*%Y
      }
      det.XhXh[i] <- det(t(X)%*%X) #determinant of Xh'Xh for ith elemental regression
      sum.abse[i] <- sum(abs(Yall - Xall%*%beta.h[,i])) #sum of absolute residuals for ith elemental regression 
    }
    callt <- match.call()
    pho <- c()
    rank.err <- rank(sum.abse)
    for (j in 1:s) { #generate pho weight for outliers
      if (rank.err[j] <= r) {
        pho[j] <- 1
      } else {
        pho[j] <- 0
      }
    }
  TEE.est <- as.matrix(t(rowSums(t(matrix(c(det.XhXh*pho,det.XhXh*pho), nrow = s, ncol = p))*beta.h)/sum(det.XhXh*pho)))
  }
  colnames(TEE.est) <- c(names)
  rownames(TEE.est) <- ""
  resid <- if (is.null(offset)) {
                 Yall - Xall%*%t(TEE.est)
                 } else {
                 Yall - (Xall%*%t(TEE.est) + offset)
                 }
  fitted <- if (is.null(offset)) {
                 Xall%*%t(TEE.est)
                 } else {
                 Xall%*%t(TEE.est) + offset
                 }
  output <- list(call = callt, formula = formula, coefficients = TEE.est, residuals = t(resid), fitted.values = t(fitted))
  return(output)
}
