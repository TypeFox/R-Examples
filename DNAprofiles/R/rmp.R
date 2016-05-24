#' Random match probability of profile(s)
#'
#' Computes the random/conditional match probability.
#' @param x Integer matrix with the profile(s) for which random match probability is computed.
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus.
#' @param markers Character vector stating the markers to use in the rmp computation. Defaults to all markers of the profile.
#' @param theta Numeric value specifying the amount of background relatedness.
#' @param cmp Logical conditional match probability. If TRUE, the Balding-Nichols formula is used to compute the conditional match probability in the same subpopulation.
#' @param ret.per.marker Logical. If TRUE, return a matrix of random match probabilities, where the columns correspond to markers.
#' @details When \eqn{\theta=0}, the simple product rule is used. Assuming Hardy-Weinberg and Linkage Equilibrium, the random match probability (rmp) for unordered is computed as the product of \deqn{2^H f_a f_b} over the loci, where \eqn{f_a} and \eqn{f_b} are respectively the population frequencies of allele \eqn{a} and \eqn{b} and \eqn{H} is the indicator function for heterozygosity (alleles \eqn{a} and \eqn{b} are not the same).
#' 
#'          When \eqn{\theta>0} and \code{cmp=FALSE}, the product rule is used that incorporates a correction for inbreeding as measured by \eqn{theta}. 
#'          
#'          When \eqn{\theta>0} and \code{cmp=TRUE}, a product rule involving a subpopulation correction is used, as given by Balding & Nichols. The match probability for homozygotes is given by: \deqn{\frac{(2 \theta+(1-\theta)f_a)(3 \theta+(1-\theta)f_a)}{(1+\theta)(1+2 \theta)},}and for heterozygotes by: \deqn{\frac{2(\theta+(1-\theta)f_a)(\theta+(1-\theta)f_b)}{(1+\theta)(1+2\theta)}.}
#'          
#'          If \code{x} contains missing values (NAs) at a marker, then the returned match probability equals one for persons with both alleles missing. If a single allele is missing, then the match probability is equal to the frequency of the single allele that is seen, unless the conditional match probability is computed.
#' @return numeric matrix of random match probabilities. When \code{ret.per.matrix} is \code{TRUE}, the columns contain rmps per marker.
#' @examples
#' ## compute the conditional match probability for two markers
#' data(freqsNLngm)
#' y <- sample.profiles(N=1,freqsNLngm)
#' rmp(y,markers = c("FGA","TH01"),theta=0.03,cmp=TRUE,ret.per.marker = TRUE)
#' rmp(y,markers = c("FGA","TH01"),ret.per.marker = TRUE) # compare to product rule estimate
#'
#' ## make a plot of density estimates of RMPs of profiles on the 10 SGMplus
#'
#' data(freqsNLsgmplus)
#' 
#' #sample profiles
#' x <- sample.profiles(N=1e3,freqsNLsgmplus)
#' 
#' #compute RMPs
#' x.rmp <- rmp(x)
#' 
#'  plot(density(log10(x.rmp)),
#'     xlab=expression(log[10](RMP)),
#'     main="Random match probabilities for SGMplus profiles")
#'     
#'
#'@export
#'
#'
rmp <- function(x,freqs=get.freqs(x),markers=get.markers(x),theta=0,cmp=FALSE,ret.per.marker=FALSE){  
  x <- Zassure.matrix(x)
  
  x.markers <- get.markers(x) # does a check on the column names of x
  freqs.markers <- names(freqs)
  
  # check if profiles and frequencies are available for these markers
  if (!all(markers %in% freqs.markers)){      
    stop("Allele frequencies unavailable for marker(s) ",paste(markers[!markers %in% freqs.markers],collapse=", "))
  }
  if (!all(markers %in% x.markers)){      
    stop("x does not contain marker(s) ",paste(markers[!markers %in% x.markers],collapse=", "))
  }
  
  if (!cmp){
    # no cmp so we can use Zrmpcpp function, which requires matrices as input
    x <- x[,paste(rep(markers,each=2),1:2,sep="."),drop=FALSE]
    # check if db contains off ladder alleles which could crash the Zrmpcpp function
    min.x <- min(x,na.rm = TRUE)
    max.x <- max(x,na.rm = TRUE)
    if (min.x<1L) stop("alleles should be positive integers")
    freqs.L <- sapply(freqs[markers],length)
    if (max.x>max(freqs.L)) stop("db contains allele that is not in freqs")
    freqs.mat <- matrix(0.,nrow = max(freqs.L),ncol = length(markers))
    for(m in seq_along(markers)) freqs.mat[seq_len(freqs.L[m]),m] <- as.vector(freqs[[markers[m]]])
    
    ret <- Zrmp(x,freqs.mat,theta,as.integer(ret.per.marker))
  }else{ # Balding-Nichols formula
    ret <- matrix(1,nrow=nrow(x),ncol=ifelse(ret.per.marker,length(markers),1))
    
    for (m in seq_along(markers)){
      ret.marker <- rep(1,nrow(x))
      
      a <- as.integer(x[,paste(markers[m],1,sep=".")]) # first allele @ marker
      b <- as.integer(x[,paste(markers[m],2,sep=".")]) # second allele @ marker      
      fr0 <- as.vector(freqs[[markers[m]]]) # look up allele freqs for marker
      
      # if a single allele is missing, then we can't use the standard formula
      a.na <- is.na(a);  b.na <- is.na(b)
      single.na <- xor(a.na,b.na)
      if (any(single.na)) {
        ret.marker[(single.na&a.na)] <- fr0[b[(single.na&a.na)]]*(1-theta)+theta
        ret.marker[(single.na&b.na)] <- fr0[a[(single.na&b.na)]]*(1-theta)+theta
      }
      
      # now the case that no allele is missing or both are missing
      hom <- (a==b) #homozygous or NA
      hom.ind <- which(hom); het.ind <- which(!hom) # gets rid of the NAs in the logical, so actual hom/hets
      
      ret.marker[hom.ind] <- (2*theta+(1-theta)*fr0[a[hom.ind]])*(3*theta+(1-theta)*fr0[a[hom.ind]])/ ((1+theta)*(1+2*theta))
      ret.marker[het.ind] <- 2*(theta+(1-theta)*fr0[a[het.ind]])*(theta+(1-theta)*fr0[b[het.ind]])/((1+theta)*(1+2*theta))
      
      if (ret.per.marker) ret[,m] <- ret.marker else ret <- ret * ret.marker
    }
  }
  
  if (ret.per.marker) colnames(ret) <- markers
  ret
}