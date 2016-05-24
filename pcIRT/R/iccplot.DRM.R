#'@rdname iccplot
#'@method iccplot DRM
#'@export

iccplot.DRM <- function(object, items="all", ...){

pp <- seq(-5,5, by=0.5)

if(is.numeric(items)){
  y <- sapply(items, function(l){

    resp <- exp(pp-object$itempar[l]) / (1+exp(pp-object$itempar[l]))
    par(ask=TRUE)
    plot(pp, resp, type="l", xlab="theta", ylab="response probability", main=paste0(names(object$itempar)[l]))
    par(ask=FALSE)
  })
} else if(items == "all"){
y <- sapply(1:length(object$itempar), function(l){

        resp <- exp(pp-object$itempar[l]) / (1+exp(pp-object$itempar[l]))
        par(ask=TRUE)
        plot(pp, resp, type="l", xlab="theta", ylab="response probability", main=paste0(names(object$itempar)[l]))
        par(ask=FALSE)
  })
}  else {stop("Items must be a numeric vector to choose a subset of items or must be 'all' to choose all items")
  }

}

