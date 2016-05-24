#' Elastic Predicition from Regression Models
#'
#' This function performs prediction from an elastic regression model
#'  with phase-variablity
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param model list describing model from elastic regression methods
#' @param y responses of test matrix f (default=NULL)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @return Returns a list containing
#' \item{y_pred}{predicited values of f or probabilities depending on model}
#' \item{SSE}{sum of squared errors if linear}
#' \item{y_labels}{labels if logistic model}
#' \item{PC}{probability of classification if logistic}
#' @keywords srvf alignment regression
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Elastic Functional Logistic Regression with Application to Physiological Signal Classification,
#'  Electronic Journal of Statistics (2014), submitted.
elastic.prediction <- function(f, time, model, y=NULL, smooth_data = FALSE, sparam = 25){
  binsize = mean(diff(time))
  eps = .Machine$double.eps

  # Compute q-function of the functional data
  tmp = gradient.spline(f,binsize,smooth_data)
  f = tmp$f
  q = tmp$g/sqrt(abs(tmp$g)+eps)

  n = ncol(q)
  method <- pmatch(model$mode, c("linear", "logistic", "mlogistic"))
  if (is.na(method))
    stop("invalid model type")
  if (method == 1 || method == 2 ){
    y_pred = rep(0,n)
  }

  if (method == 3){
    m = model$n_classes
    y_pred = matrix(0,n,m)
  }

  for (ii in 1:n){
    diff1 = model$q - replicate(ncol(model$q), q[,ii])
    dist = apply(abs(diff1)^2,2,sum)^(1/2)
    argmin = which.min(dist)
    q_tmp = warp_q_gamma(time, q[,ii], model$gamma[, argmin])
    if (method == 1){
      y_pred[ii] = model$alpha + trapz(time, q_tmp*model$beta)
    } else if (method ==2){
      y_pred[ii] = model$alpha + trapz(time, q_tmp*model$beta)
    } else if (method == 3){
      for (jj in 1:m){
        y_pred[ii,jj] = model$alpha[jj] + trapz(time, q_tmp*model$beta[,jj])
      }
    }
  }

  if (is.null(y)){
    if (method == 1){
      SSE = NULL
    } else if (method == 2){
      y_pred = phi(y_pred)
      y_labels = rep(1,n)
      y_labels[y_pred<0.5] = -1
      PC = NULL
    } else if (method == 3){
      y_pred = phi(as.vector(y_pred))
      y_pred = array(y_pred,c(n, m))
      y_labels = apply(y_pred, 1, which.max)
      PC = NULL
    }
  } else {
    if (method == 1){
      SSE = sum((y-y_pred)^2)
    } else if (method == 2){
      y_pred = phi(y_pred)
      y_labels = rep(1,n)
      y_labels[y_pred<0.5] = -1
      TP = sum(y[y_labels == 1] == 1)
      FP = sum(y[y_labels == -1] == 1)
      TN = sum(y[y_labels == -1] == -1)
      FN = sum(y[y_labels == 1] == -1)
      PC = (TP+TN)/(TP+FP+FN+TN)
    } else if (method == 3){
      y_pred = phi(as.vector(y_pred))
      y_pred = array(y_pred,c(n, m))
      y_labels = apply(y_pred, 1, which.max)
      PC = rep(0,m)
      cls_set = 1:m
      for (ii in 1:m){
        cls_sub = setdiff(cls_set,ii);
        TP = sum(y[y_labels == ii] == ii)
        FP = sum(y[is.element(y_labels, cls_sub)] == ii)
        TN = sum(y[is.element(y_labels, cls_sub)] ==
                   y_labels[is.element(y_labels, cls_sub)])
        FN = sum(is.element(y[y_labels == ii], cls_sub))
        PC[ii] = (TP+TN)/(TP+FP+FN+TN)
      }
      PC = sum(y == y_labels)/length(y_labels)
    }
  }

  if (method == 1){
    out = list(y_pred=y_pred,SSE=SSE)
  } else if (method == 2){
    out = list(y_prob=y_pred, y_labels=y_labels, PC=PC)
  } else if (method == 3){
    out = list(y_prob=y_pred, y_labels=y_labels, PC=PC)
  }

  return(out)

}
