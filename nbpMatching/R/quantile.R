#'Quantile for upper-triangular values in distance matrix
#'
#'Extend the stats \code{\link[stats]{quantile}} function for handling distancematrix objects.
#'
#'The upper.triangular values of the distance matrix object are passed to the
#'\code{\link{quantile}} function.
#'
#'@aliases quantile quantile,distancematrix-method
#'@param x A distancematrix object.
#'@param probs numeric vector or probabilities with values in [0,1].
#'@param \dots Additional arguments, passed to \code{\link[stats]{quantile}}.
#'@return numeric vector of quantiles corresponding to the given probabilities
#'@export
#'@rdname quantile
#'@author Cole Beck
#'@examples
#'
#'plainmatrix<-as.matrix(dist(sample(1:25, 8, replace=TRUE)))
#'mdm<-distancematrix(plainmatrix)
#'quantile(mdm, probs=c(0.0, 0.25, 0.50, 0.75, 1.00))
#'

setMethod("quantile", "distancematrix", function(x, probs, ...) {
    if(missing(probs)) probs=seq(from=0, to=100)/100
    quantile(x[upper.tri(x)], probs=probs, ...)
})
