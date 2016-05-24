W.combine <- function(W,n) {
  if(!is.list(W) || !is.numeric(n))
    stop("W must be a list of W values to combine, with n a vector of sample sizes.\n")
  if(length(W) != length(n))
    stop("need equal length W and n.\n")
  for(i in 1:length(W))
    W[[i]] <- W[[i]]/(n[[i]] + 1)
  if(length(W[[1]])>1)
    return(rowSums(do.call("cbind",W)))
  ##    W <- mapply(function(x,y) { x/y }, W, as.list(n))
  ##    W <- unlist(W)/n
  return(sum(unlist(W)))
}




##' Calculate a Z score from a Wilcoxon statistic and a set of random Wilcoxon
##' statistics
##' 
##' The mean of a Wilcoxon statistic is unaffected by correlation within the
##' variable under test, but its variance is.  This function uses a set of
##' Wilcoxon statistics generated from permuted data to estimate the variance
##' empirically, and thus calculate a Z score.
##' 
##' 
##' @param W Wilcoxon statistic for observed data.
##' @param Wstar A vector of Wilcoxon statistics for a set of permuted data.
##' @param n.in The number of items (SNPs) in the regions to be tested.
##' @param n.out The number of items (SNPs) in the control regions.
##' @return A list with two elements:
##' \describe{
##'
##' \item{Z.theoretical}{which uses the theoretical mean of the Wilcoxon
##' distribution under the null generated from n.in, n.out above}
##'
##' \item{Z.empirical}{which uses Wstar to calculate an empirical estimate
##' of the mean of the Wilcoxon distribution under the null}
##'
##' }
##' @note The function can also deal with combining W statistics from multiple
##' strata, as is typical in a meta analysis of GWAS data, using van Elteren's
##' method.  Strata may be defined by different geography or different SNP
##' chips.
##' @author Chris Wallace
##' @export
##' @seealso \code{\link{wilcoxon}}
##' @keywords htest
##' @examples
##' 
##' x <- exp(-rexp(1000)) # uniform
##' y <- exp(-rexp(1000,0.8)) # skewed towards 0
##' W <- wilcoxon(p=c(x,y),snps.in=1:1000)
##' 
##' p.perm <- matrix(sample(c(x,y),replace=TRUE,size=10000),ncol=5)
##' Wstar <- wilcoxon(p=p.perm,snps.in=1:1000)
##' 
##' Z.value(W=W, Wstar=Wstar, n.in=1000, n.out=1000)
##' 
Z.value <- function(W,Wstar,n.in,n.out) {
  if(is.list(W) && (!is.list(Wstar) || length(W)!=length(Wstar) || length(W)!=length(n.in) || length(W) !=length(n.out)))
    stop("W, Wstar must each be lists of length equal to the length of n.in and n.out.\n")
  if(is.list(W)) {
    W <- W.combine(W,n.in+n.out)
    Wstar <- W.combine(Wstar,n.in+n.out)
    mu.theor <- sum(mu.wilcox(n.in,n.out) / (n.in+n.out+1))
  } else {
    mu.theor <- mu.wilcox(n.in,n.out)
  }
  sd.est <- sd(Wstar)
  mu.est <- mean(Wstar)
  z <- (W-mu.theor)/sd.est
  Z.theoretical <- list(statistic=c(Z=z),
                        p.value=2*pnorm(abs(z),lower.tail=FALSE),
#                        parameter=c(mean.null=mu.theor),
                        method="Wilcoxon theoretical mean",
                        data.name=deparse(substitute(W)))
  z <- (W-mu.est)/sd.est
  Z.empirical <- list(statistic=c(Z=z),
                      p.value=2*pnorm(abs(z),lower.tail=FALSE),
#                      parameter=c(mean.null=mu.est),
                      method="Wilcoxon empirical mean",
                      data.name=deparse(substitute(W)))
  class(Z.theoretical) <- "htest"
  class(Z.empirical) <- "htest"
  return(list(Z.theoretical=Z.theoretical,
              Z.empirical=Z.empirical))
}
