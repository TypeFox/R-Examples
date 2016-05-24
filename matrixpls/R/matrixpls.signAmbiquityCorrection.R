# =========== Weight sign corrections ===========

#'Sign ambiquity corrections
#'
#'Sign ambiquity corrections adjust the signs of the weights to satisfy a criterion.
#'
#'Instead of fixing a weight to a particular value, composite variables are typically provided a
#'scale by standardization. This leads to sign indeterminacy because standardized weights \code{W}
#'and \code{-W} both satisfy the scalign constraint. The sing ambiquity corrections add additional
#'constraints that make 
#'
#'The sign indeterminacy
#'corrections should not be confused with sign chance corrections applied to boostrap samples 
#'(See \code{\link{signChangeCorrection}}).
#'
#'@inheritParams matrixpls-common
#'@return \code{W} after sign correction.
#'
#'@name signAmbiquityCorrection
#'
#'@references
#'Wold, H. (1985). Partial Least Squares. In S. Kotz & N. L. Johnson (Eds.), Encyclopedia  of 
#'statistical sciences (Vol. 6, pp. 581–591). New York: Wiley.

#'@seealso 
#'\code{\link{matrixpls}};
#
NULL

#'@describeIn signAmbiquityCorrection Adjust the signs of W so that the majority of the indicators are positively 
#'correlated with the composite as proposed by Wold (1985).
#'
#'@export

weightSign.Wold1985 <- function(W,S){
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  signSums <- rowSums(sign(IC * (W!=0)))
  sweep(W, 1,ifelse(signSums<0, -1, 1),"*")
}

#'@describeIn signAmbiquityCorrection Adjust the signs of W so that the first indicator of each composite has positive
#'weight.
#'
#'@export

weightSign.dominantIndicator <- function(W,S){
  # Calculate the covariance matrix between indicators and composites
  signs <- apply(W,1,function(x){
    i <- min(which(x!=0))
    sign(x[i])
  })
  sweep(W, 1,signs,"*")
}
