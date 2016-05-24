#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: kldiv.R 4834 2012-08-02 10:17:09Z gruen $
#

setMethod("KLdiv", "matrix",
function(object, eps=10^-4, overlap=TRUE,...)
{
    if(!is.numeric(object))
        stop("object must be a numeric matrix\n")
    
    z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
    colnames(z) <- rownames(z) <- colnames(object)
    
    w <- object < eps
    if (any(w)) object[w] <- eps
    object <- sweep(object, 2, colSums(object) , "/")
    
    for(k in seq_len(ncol(object)-1)){
      for(l in 2:ncol(object)){
        ok <- (object[, k] > eps) & (object[, l] > eps)
        if (!overlap | any(ok)) {
          z[k,l] <- sum(object[,k] *
                        (log(object[,k]) - log(object[,l])))
          z[l,k] <- sum(object[,l] *
                        (log(object[,l]) - log(object[,k])))
        }
      }
    }
    diag(z)<-0
    z
})

setMethod("KLdiv", "flexmix",
function(object, method = c("continuous", "discrete"), ...) {
  method <- match.arg(method) 
  if (method == "discrete") z <- KLdiv(object@posterior$scaled, ...)
  else {
    z <- matrix(0, object@k, object@k)
    for (i in seq_along(object@model)) {
      comp <- lapply(object@components, "[[", i)
      z <- z + KLdiv(object@model[[i]], comp)
    }
  }
  z
})

setMethod("KLdiv", "FLXMRglm",
function(object, components, ...) {
  z <- matrix(NA, length(components), length(components))
  mu <- lapply(components, function(x) x@predict(object@x))
  if (object@family == "gaussian") {
    sigma <- lapply(components, function(x) x@parameters$sigma)
    for (k in seq_len(ncol(z)-1)) {
      for (l in seq_len(ncol(z))[-1]) {
        z[k,l] <- sum(log(sigma[[l]]) - log(sigma[[k]]) + 1/2 * (-1 + ((sigma[[k]]^2 + (mu[[k]] - mu[[l]])^2))/sigma[[l]]^2))
        z[l,k] <- sum(log(sigma[[k]]) - log(sigma[[l]]) + 1/2 * (-1 + ((sigma[[l]]^2 + (mu[[l]] - mu[[k]])^2))/sigma[[k]]^2))
      }
    }
  }
  else if (object@family == "binomial") {
    for (k in seq_len(ncol(z)-1)) {
      for (l in seq_len(ncol(z))[-1]) {
        z[k,l] <- sum(mu[[k]] * log(mu[[k]]/mu[[l]]) + (1-mu[[k]]) * log((1-mu[[k]])/(1-mu[[l]])))
        z[l,k] <- sum(mu[[l]] * log(mu[[l]]/mu[[k]]) + (1-mu[[l]]) * log((1-mu[[l]])/(1-mu[[k]])))
      }
    }
  }
  else if (object@family == "poisson") {
    for (k in seq_len(ncol(z)-1)) {
      for (l in seq_len(ncol(z))[-1]) {
        z[k,l] <- sum(mu[[k]] * log(mu[[k]]/mu[[l]]) + mu[[l]] - mu[[k]])
        z[l,k] <- sum(mu[[l]] * log(mu[[l]]/mu[[k]]) + mu[[k]] - mu[[l]])
      }
    }
  }
  else if (object@family == "gamma") {
    shape <- lapply(components, function(x) x@parameters$shape)
    for (k in seq_len(ncol(z)-1)) {
      for (l in seq_len(ncol(z))[-1]) {
        X <- mu[[k]]*shape[[l]]/mu[[l]]/shape[[k]]
        z[k,l] <- sum(log(gamma(shape[[l]])/gamma(shape[[k]])) + shape[[l]] * log(X) - shape[[k]] * (1 - 1/X) +
                      (shape[[k]] - shape[[l]])*digamma(shape[[k]]))
        z[l,k] <- sum(log(gamma(shape[[k]])/gamma(shape[[l]])) - shape[[k]] * log(X) - shape[[l]] * (1 - X) +
                      (shape[[l]] - shape[[k]])*digamma(shape[[l]]))
      }
    }
  }
  else stop(paste("Unknown family", object@family))
  diag(z) <- 0
  z
})

setMethod("KLdiv", "FLXMC",
function(object, components, ...) {
  z <- matrix(NA, length(components), length(components))
  if (object@dist == "mvnorm") {
    center <- lapply(components, function(x) x@parameters$center)
    cov <- lapply(components, function(x) x@parameters$cov)
    for (k in seq_len(ncol(z)-1)) {
      for (l in seq_len(ncol(z))[-1]) {
        z[k,l] <- 1/2 * (log(det(cov[[l]])) - log(det(cov[[k]])) - length(center[[k]]) +
                         sum(diag(solve(cov[[l]]) %*% (cov[[k]] + tcrossprod(center[[k]] - center[[l]])))))
        z[l,k] <- 1/2 * (log(det(cov[[k]])) - log(det(cov[[l]])) - length(center[[l]]) +
                         sum(diag(solve(cov[[k]]) %*% (cov[[l]] + tcrossprod(center[[l]] - center[[k]])))))
      }
    }
  }
  else stop(paste("Unknown distribution", object@dist))
  diag(z) <- 0
  z
})

###**********************************************************


