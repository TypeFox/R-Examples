#'Sign change corrections for bootstrap
#'
#'These functions selectively reverse the signs of the weights in boostrap samples to be consistent 
#'with the weights calculated based on the original sample.
#'
#'Sign change corrections are a controversial and inconsistently implemented feature in PLS analysis.
#'The two corrections described in the literature are the individual sign chance correction and the
#'construct level sign chance corrections.

#'The individual correction changes the signs of \code{W} to match \code{Worig}.
#'
#'The construct level correction changes the signs of \code{W} on all rows where the sign of the
#'sum of the row differs between \code{Worig} and \code{W}.
#'
#'The sign chance corrections are described ambiquosly and sometimes implemented inconsistently
#'between software. \pkg{matrixpls} implements the corrections by adjusting the weights before
#'calculating parameter estimates in each bootstrap replication. Some software implement the
#'correction post-hoc by adjusting the bootstrap estimates directly. Moreover, the 
#'literature present at least two different formulas for the construct level correction. 
#'\pkg{matrixpls} implements the version described by Tenenhaus et al. (2005).
#'
#'The sign chance
#'corrections should not be confused with sign indeterminacy corrections applied to 
#'individual analyses
#'(See \code{\link{signAmbiquityCorrection}}).

#'
#'@param Worig The original weight matrix.
#'
#'@param W a Weight matrix of a bootstrap sample.
#'
#'@return A weight matrix with the same dimensions as \code{W} after applying the correction.
#'
#'@seealso
#'\code{\link{matrixpls.boot}}
#'
#'@name signChangeCorrection
#'
#'
#'@references 
#'
#'Tenenhaus, M., Esposito Vinzi, V., Chatelin, Y.-M., & Lauro, C. (2005). PLS Path Modeling.
#'\emph{Computational Statistics & Data Analysis}, 48(1), 159–205. doi:10.1016/j.csda.2004.03.005
#'
#'Rönkkö, M., McIntosh, C. N., & Antonakis, J. (2015). On the adoption of partial least squares in 
#'psychological research: Caveat emptor. \emph{Personality and Individual Differences}, (87), 76–84.
#'\href{https://doi.org/10.1016/j.paid.2015.07.019}{DOI:10.1016/j.paid.2015.07.019}
NULL

#'@describeIn signChangeCorrection individual sign change correction
#'@export


signChange.individual <- function(Worig,W){
  W * sign(Worig) * sign(W)
}

#'@describeIn signChangeCorrection individual sign change correction
#'@export

signChange.construct <- function(Worig,W){
  origSums <- apply(Worig,1,sum)
  curSums <- apply(W,1,sum)
  signs <- ifelse(abs(origSums-curSums) > abs(origSums+curSums), -1,1)
  sweep(W,MARGIN=1,signs,`*`)
}
