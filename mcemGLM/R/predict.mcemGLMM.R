predict.mcemGLMM <- function(object, newdata, type = c("link", "response"), se.fit = FALSE, ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  if (missing(newdata)) {
    lin0 <- as.vector(object$x %*% coef0)
    if (type[1] == "link") {
      if (se.fit == FALSE)
        return(lin0)
      else {
        kN <- length(object$y)
        tmp <- matrix(0, kN, 2)
        tmp[, 1] <- lin0
        colnames(tmp) <- c("Estimate", "SE")
        cmat <- covMat.mcemGLMM(object)
        for (i in 1:kN) {
          tmp[i, 2] <- object$x[i, ] %*% cmat %*% (object$x[i, ])
        }
        return(tmp)
      }
    } 
    if (type == "response") {
      if (object$call$family == "bernoulli") {
        if (se.fit == FALSE)
          return(exp(lin0) / (1 + exp(lin0)))
        else {
          kN <- length(object$y)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0) / (1 + exp(lin0))
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- covMat.mcemGLMM(object)
          for (i in 1:kN) {
            tmp[i, 2] <- exp(2 * lin0[i])/(1 + exp(lin0[i]))^4 * object$x[i, ] %*% cmat %*% object$x[i, ]
          }
          return(tmp)
        }
      }
      if (object$call$family %in% c("poisson", "negbinom", "gamma")) {
        if (se.fit == FALSE)
          return(exp(lin0))
        else {
          kN <- length(object$y)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0)
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- covMat.mcemGLMM(object)
          for (i in 1:kN) {
            tmp[i, 2] <- exp(2 * lin0[i]) * object$x[i, ] %*% cmat %*% object$x[i, ]
          }
          return(tmp)
        }
      }
    }
  } else {
    # newdata
    tmp.x <- model.matrix(as.formula(object$call$fixed)[-2], data = newdata)
    if (!prod(colnames(tmp.x) == colnames(object$x))) {
      stop("Incorrect new data.")
    }
    lin0 <- as.vector(tmp.x %*% coef0)
    if (type[1] == "link") {
      if (se.fit == FALSE)
        return(lin0)
      else {
        kN <- nrow(tmp.x)
        tmp <- matrix(0, kN, 2)
        tmp[, 1] <- lin0
        colnames(tmp) <- c("Estimate", "SE")
        cmat <- covMat.mcemGLMM(object)
        for (i in 1:kN) {
          tmp[i, 2] <- tmp.x[i, ] %*% cmat %*% (tmp.x[i, ])
        }
        return(tmp)
      }
    } 
    if (type == "response") {
      if (object$call$family == "bernoulli") {
        if (se.fit == FALSE)
          return(exp(lin0) / (1 + exp(lin0)))
        else {
          kN <- nrow(tmp.x)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0) / (1 + exp(lin0))
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- covMat.mcemGLMM(object)
          for (i in 1:kN) {
            tmp[i, 2] <- exp(2 * lin0[i])/(1 + exp(lin0[i]))^4 * tmp.x[i, ] %*% cmat %*% tmp.x[i, ]
          }
          return(tmp)
        }
      }
      if (object$call$family %in% c("poisson", "negbinom", "gamma")) {
        if (se.fit == FALSE)
          return(exp(lin0))
        else {
          kN <- nrow(tmp.x)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0)
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- covMat.mcemGLMM(object)
          for (i in 1:kN) {
            tmp[i, 2] <- exp(2 * lin0[i]) * tmp.x[i, ] %*% cmat %*% tmp.x[i, ]
          }
          return(tmp)
        }
      }
    }
  }
}