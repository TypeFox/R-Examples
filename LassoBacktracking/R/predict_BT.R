#' Make predictions from a "\code{BT}" object.
#' 
#' Similar to other predict methods, this function predicts fitted values and computes coefficients
#' from a fitted "\code{BT}" object.
#' 
#' @param object fitted "\code{BT}" object.
#' @param newx matrix of new values of design matrix at which predictions are to be made. Ignored
#'  when \code{type=="coefficients"}.
#' @param s value of the penalty parameter at which predictions are required. If the value
#'  is not one of the \code{lambda} values present in \code{object} the output will be
#'  etermined by linear interpolation. Default is the entire sequence of \code{lambda}
#'  values present in \code{object}.
#' @param iter iteration at which predictions are required. Default is the entire sequence
#'  of iterations in \code{object}.
#' @param type of prediction required. Type "\code{response}" gives estimates of the response
#'  whilst type "\code{coefficients}" gives coefficient estimates.
#' @param ... not used. Other arguments to \code{predict}.
#' @return Either a vector of predictions or, if either \code{s} or \code{iter} are \code{NULL},
#'  a three-dimensional array with last two dimensions indexing different \code{lambda} values and
#'  iterations.
#' @examples
#' x <- matrix(rnorm(100*250), 100, 250)
#' y <- x[, 1] + x[, 2] - x[, 1]*x[, 2] + x[, 3] + rnorm(100)
#' out <- LassoBT(x, y, iter_max=10)
#' predict(out, newx=x[1:2, ])
#' @export
predict.BT <- function(object, newx, s=NULL, iter=NULL, type=c("response", "coefficients"), ...) {
  max_iter <- ncol(object$path_lookup)
  nlambda <- length(object$lambda)
  p <- object$nvars
  
  if (!missing(newx)) {
    if (!is.matrix(newx)) stop("newx must be a matrix")
    if (ncol(newx) < object$var_indices[p])
      stop(paste("newx must have at least", object$var_indices[p], "variables"))
    n <- nrow(newx)
  } else {
    n <- object$nobs
  }
  
  type <- match.arg(type)
  
  # if exactly one of s or iter is not given, we ignore the other that is given
  if (is.null(s) || is.null(iter)) {
    if (type != "response") stop("if either s or iter are NULL, type must be response")
    #dimnames=c("observations", "lambda", "iterations")
    out <- array(dim=c(n, nlambda, max_iter))
    n_inter <- ncol(object$interactions)
    if (!missing(newx)) {
      x_all <- matrix(nrow=n, ncol=p + n_inter)
      k <- p
      for (j in seq_len(n_inter)) {
        k <- k + 1L
        x_all[, k] <- newx[, object$interactions[1L, j]]*newx[, object$interactions[2L, j]]
      }
      x_all[, 1L:p] <- newx[, object$var_indices]
      object$fitted <- lapply(seq_along(object$beta),
                              function(k) as.matrix(x_all %*% object$beta[[k]] + rep(object$a0[[k]], each=n)))
    }
    # could code this in c++ for speed
    for (iter in seq_len(max_iter)) {
      for (l in seq_len(nlambda)) {
        path <- object$path_lookup[l, iter]
        out[, l, iter] <- object$fitted[[path]][, l - object$l_start[path] + 1L]
      }
    }
    return(out)
  }

  iter <- as.integer(iter)
  if (length(iter) > 1L) stop("iter must be a single number")
  if (iter < 1L) {
    iter <- 1L
    warning("iter set to 1")
  } else if (iter > max_iter) {
    warning(paste("iter outside range so set to", max_iter))
    iter <- max_iter        
  }
  
  s <- as.double(s)
  if (length(s) > 1L) stop("s must be a single number")
  out_interp <- interp_w(s, object$lambda)
  l1 <- out_interp$indices[1L]
  l2 <- out_interp$indices[2L]
  w1 <- out_interp$weights[1L]
  w2 <- out_interp$weights[2L]
  path1 <- object$path_lookup[l1, iter]
  path2 <- object$path_lookup[l2, iter] 
  l1 <- l1 - object$l_start[path1] + 1L
  l2 <- l2 - object$l_start[path2] + 1L
  beta1 <- object$beta[[path1]][, l1]
  beta2 <- object$beta[[path2]][, l2]
  mu1 <- object$a0[[path1]][l1]
  mu2 <- object$a0[[path2]][l2]
  
  if (type == "coefficients") {
    warning("newx is not needed to calculate coefficients")
    coef_out <- c(w1*mu1 + w2*mu2, w1*beta1 + w2*beta2)
    names(coef_out) <- c("intercept", rownames(object$beta[[path1]]))
    return(coef_out)
  }
  
  # type must equal response
  if (missing(newx)) {
    if (type == "response") {
      return(as.numeric(w1*object$fitted[[path1]][, l1] +
                          w2*object$fitted[[path2]][, l2]))
    }
  }
  
  if (type == "response") {
    nonzero_inter1 <- which(beta1 != 0)
    nonzero_inter1 <- nonzero_inter1[nonzero_inter1 > p]
    nonzero_inter2 <- which(beta2 != 0)
    nonzero_inter2 <- nonzero_inter2[nonzero_inter2 > p]
    nonzero_inter <- unique(c(nonzero_inter1, nonzero_inter2))
    
    if (length(nonzero_inter) > 0L) {
      x_inter <- matrix(nrow=n, ncol=length(nonzero_inter))
      for (j in seq_along(nonzero_inter)) {
        x_inter[, j] <- newx[, object$interactions[1, nonzero_inter[j]-p]]*
          newx[, object$interactions[2, nonzero_inter[j]-p]]
      }
      inter_sig <- as.numeric(x_inter %*% (w1*beta1[nonzero_inter] + w2*beta2[nonzero_inter]))
    } else {
      inter_sig <- 0
    }
    #x_all[, 1L:p] <- newx[, object$var_indices]
    #var_all <- c(1:object$nvars, nonzero_inter)
    return(as.numeric(newx[, object$var_indices] %*% (w1*beta1[1:p] + w2*beta2[1:p])) +
           inter_sig + w1*object$a0[[path1]][l1] + w2*object$a0[[path2]][l2])
  }
}

#' @rdname predict.BT
#' @export
coef.BT <- function(object, s=NULL, iter=NULL, ...) {
  return(predict.BT(object = object, s=s, iter=iter, ...))
}