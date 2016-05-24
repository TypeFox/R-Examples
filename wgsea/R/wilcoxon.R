mu.wilcox <-
function(n1,n2) {
  exp(log(n1)+log(n2)-log(2))
}
sd.wilcox <- function(n1,n2) {
  sqrt( n1*(n2) * (n1+n2+1)/12 )
}

calc.wilcoxon <-
function(p,snps.in,n0,w=NULL) {
  R <- rank(p)[snps.in]
  if(!is.null(w))
    R <- R * w
  return(sum(R) - n0 * (n0+1)/2)
}





##' Wilcoxon test statistic, with optional weights.
##' 
##' Calculate a Wilcoxon two group rank test statistic, with optional
##' propensity score weighting.
##' 
##' 
##' @param p a numeric vector of observed p values from a list of SNPs or a
##' matrix, with each column representing a vector under a different
##' permutation of the dataset.  
##' @param snps.in a numeric vector indicating which members of p form the test
##' group (their complement form the control group). 
##' @param weights optional propensity score weights.  These are binned
##' according to binsize, and a weight calculated which is inversely
##' proportional to the probability of sampling a member of the test group in
##' that bin. 
##' @param binsize see weights, above.
##' @return A numeric value or, if p is a matrix, a numeric vector.
##' @author Chris Wallace
##' @export
##' @seealso \code{\link{Z.value}}
##' @references Propensity weights are described
##' 
##' Rosenbaum, P. R. & Rubin, D. B. The central role of the propensity score in
##' observational studies for causal effects. Biometrika, 1983, 70, 41-55
##' 
##' Rosenbaum, P. R. Model-based direct adjustment. Journal of the American
##' Statistical Association, 1987, 82, 387-394
##' @keywords htest
##' @examples
##' 
##' x <- exp(-rexp(1000)) # uniform
##' y <- exp(-rexp(1000,0.8)) # skewed towards 0
##' wilcoxon(p=c(x,y),snps.in=1:1000)
##' 
##' ## note, should be equal to
##' wilcox.test(x,y)
##' 
wilcoxon <-
function(p,snps.in,weights=NULL,binsize=0.05) {
  n0 <- length(snps.in)
  if(!is.null(weights)) {
    n <- length(weights)
    seq.cut<-seq(min(weights),max(weights),by=binsize)
    if(max(seq.cut)<=max(weights))
    	seq.cut[ which.max(seq.cut) ] <- max(weights) + binsize/100
    	
    if(min(seq.cut)>=min(weights))
    	seq.cut[ which.min(seq.cut) ] <- min(weights) - binsize/100 
    bin <- cut(weights,seq.cut)
    t.bin<-table(bin)
    if(min(t.bin)<5)
    	stop("Data is too sparse for the weights and binsize parameters input")
    f <- t.bin/n
    f0 <- table(bin[snps.in])/n0
    bin.in <- bin[snps.in]
    w <- f[bin.in]/f0[bin.in]
    if(sum(is.na(w))){
      ## unlikely that we reach here but in case we do and to prevent an
      ## score of NA being returned silently throw an error.
    	stop("Propensity score and selected binsize incompatible; NA's generated")
    }
  } else {
    w <- NULL
  }
  
  if(is.matrix(p)) { # assume permutations in each column
    nc <- ncol(p)
    W <- numeric(nc)
    for(j in 1:nc) {
      #R <- rank(p[,j])[snps.in]
      #if(!is.null(weights))
      #  R <- R * w
      W[j] <- calc.wilcoxon(p[,j],snps.in,n0,w=w)
    }
  } else {
    W <- calc.wilcoxon(p,snps.in,n0,w=w)
  }
  return(W)
}

