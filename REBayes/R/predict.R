#' Predict Method for GLmix
#' 
#' Predict Method for Gaussian Location Mixtures
#' 
#' The predict method for \code{GLmix} objects will compute means, medians or
#' modes of the posterior according to the \code{Loss} argument.  Typically,
#' \code{newdata} would be passed to \code{predict}
#' 
#' @param object fitted object of class "GLmix"
#' @param newdata Values at which prediction is desired
#' @param Loss Loss function used to generate prediction
#' @param ... optional arguments to predict
#' @return A vector of predictions
#' @author Roger Koenker
#' @keywords nonparametric
#' @export
predict.GLmix <- function(object, newdata, Loss = 2, ...) {
    x <- newdata
    v <- object$x
    fv <- object$y

    if(Loss == 2) { # mean case equivalent to object$dy when x == original data
	A <- dnorm(outer(x, v, "-"), sd = object$sigma)
	xhat <- as.vector((A %*% (fv * v))/(A %*% fv))
    }
    else if(Loss == 1){ #median case
       A <- t(t(dnorm(outer(x, v, "-"), sd = object$sigma)) * fv)
       B <- apply(A/apply(A,1,sum),1,cumsum) < 1/2
       j <- apply(B,2,sum)
       if(any(j == 0)) { # Should only happen when v grid is very restricted
	   j <- j + 1
	   warning("zeros in posterior median indices")
       }
       xhat <- v[j]
    }
    else if(Loss == 0) { # mode case
       A <- t(t(dnorm(outer(x, v, "-"), sd = object$sigma)) * fv)
       xhat <- v[apply(A/apply(A,1,sum),1,which.max)]
    }
    else 
	stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}
