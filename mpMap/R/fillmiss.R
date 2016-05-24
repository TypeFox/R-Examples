#' Fill in missing values for an mpcross object
#' 
#' Use multi-point founder probabilities to fill in the most likely value for a missing genotype
#' @export
#' @param object Object of class \code{mpcross}
#' @param threshold Threshold for probability to call an allele
#' @param dart Flag for whether genotypes are coded as 0/1; shortcut for computation if true
#' @return An mpcross object with a new component $missfinals for the original set of genotypes. The component $finals is replaced by the imputed values. 
#' @seealso \code{\link[mpMap]{mpprob}}, \code{\link[mpMap]{mpcross}}

fillmiss <- function(object, threshold=.7, dart=FALSE)
{
  if (!inherits(object, "mpprob")) stop("Must compute founder probabilities first")

  chr <- names(object$map)

  n.founders <- nrow(object$founders)
  output <- object
  output$prob <- NULL
  output$estfnd <- NULL

  probs <- do.call("cbind", object$prob)
  ifmat <- matrix(nrow=nrow(probs), ncol=ncol(probs)/n.founders)
  founder <- as.vector(object$founders)
 
  missfx <- function(x) {
	tmp <- by(x, fd, sum)
	if (max(tmp)>threshold) return(as.numeric(names(which.max(tmp)))) else return(NA)}

  for (jj in seq(1, ncol(probs), n.founders))
  {
     index <- (jj-1)/n.founders+1
     fd <- founder[jj:(jj+n.founders-1)]

     if (!dart){
      alltab <- apply(probs[,jj:(jj+n.founders-1)], 1, missfx)
      ifmat[,index] <- alltab} else {
	tmp <- probs[,jj:(jj+n.founders-1)]%*%founder[jj:(jj+n.founders-1)]
	ifmat[,index] <- as.numeric(tmp>threshold)
	}
  }


  rownames(ifmat) <- rownames(object$finals)
  colnames(ifmat) <- colnames(object$finals)
  output$missfinals <- output$finals
  output$finals <- ifmat
  attr(output$finals, "threshold") <- threshold

  # probably need to find a way to make this faster. for 0's and 1's just take
  # matrix product of probs and founders, but this doesn't work for other alleles. 
  # o.w can expect about 2000 s. 

  return(output)
}
