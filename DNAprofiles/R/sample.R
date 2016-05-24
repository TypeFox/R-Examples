#'Sample random unrelated profiles
#'
#' @param N number of profiles to sample (integer).
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param markers A character vector naming the markers of the resulting sample. Default to all markers of the \code{freqs} argument.
#' @param theta numeric value specifying the amount of background relatedness, i.e. the probability that both alleles at a locus are identical by descent.
#' @details The function randomly samples DNA profiles according to the supplied allelic frequencies.
#' 
#'          When \eqn{\theta=0}, the function assumes HW-equilibrium, so the alleles of a person at a locus are independent samples.
#'          
#'          When \eqn{\theta>0}, the alleles of a person at a locus are ibd with probability \eqn{\theta}.
#' @return An object of class \code{profiles}, which is an integer matrix with \eqn{N} rows and twice the number of loci columns. The integers correspond to the index in the allelic frequency vector, NOT to STRs. Each row is a profile and every two columns contain the two alleles at a locus.
#' @seealso \code{\link{sample.pairs}}, \code{\link{sample.relatives}}
#' @examples
#' data(freqsNLsgmplus)
#' db <- sample.profiles(N=1e3,freqs=freqsNLsgmplus)
#' @export
sample.profiles <- function(N, freqs,markers=names(freqs),theta=0){  
  # checks
  if (N<0) stop("N should not be negative")
  Zchecktheta(theta)
  
  freqs.markers <- names(freqs)
  # check if profiles and frequencies are available for these markers
  if (!all(markers %in% freqs.markers)){      
    stop("Allele frequencies unavailable for marker(s) ",paste(markers[!markers %in% freqs.markers],collapse=", "))
  }

  ret <-  matrix(integer(),nrow=N,ncol=2*length(markers))
  colnames(ret) <- paste(rep(markers,each=2),1:2,sep=".")
  
  for (m in markers){
    fr <- as.vector(freqs[[m]]); fr.n <- length(fr)
    
    ret[,paste(m,1,sep=".")] <- sample.int(fr.n, size=N, replace=TRUE, fr)
    ret[,paste(m,2,sep=".")] <- sample.int(fr.n, size=N, replace=TRUE, fr)
    
    # ibd alleles
    if (theta>0){
      ibd <- which(sample(c(TRUE,FALSE),size=N,replace=TRUE,c(theta,1-theta)))
      ret[ibd,paste(m,2,sep=".")] <- ret[ibd,paste(m,1,sep=".")]
    }    
  }
  
  class(ret) <- c("profiles",class(ret))
  attr(ret,"freqs") <- freqs
  ret 
}
NULL
#' Sample random relatives of one profile or many profiles
#' 
#' @param x An integer matrix specifying a single profile. Alternatively an integer vector containing a single profile, e.g. obtained when a row is selected from a matrix of profiles.
#' @param N number of relatives to sample per profile (integer).
#' @param type A character string giving the type of relative. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param markers A character vector naming the markers of the resulting sample. Default to all markers of the \code{freqs} argument.
#' @param theta numeric value specifying the amount of background relatedness.
#' @details When \code{x} is a single profile, the function samples \eqn{N} profile that are related to \code{x} with the supplied type of relationship (\code{type}).
#' 
#' When \code{x} is a database of profiles, there are \eqn{N} relatives sampled per profile in \eqn{x}. Hence there will be \eqn{nrow(x)*N} profiles returned. In this case the returned matrix contains the relatives in order of the profiles they correspond with, e.g. the relatives of 1,1,2,2,3,3,.. when \eqn{N=2}.
#' @return An object of class \code{profiles}, which is an integer matrix with \eqn{N*nrow(x)} rows and twice the number of loci columns. The integers correspond to the index in the allelic frequency vector, NOT to STRs. Each row is a profile and every two columns contain the two alleles at a locus.
#' @seealso \code{\link{sample.pairs}}, \code{\link{sample.profiles}}
#' @examples
#' ## sample either many relatives of one profile, or one (or more) relatives of many profiles
#' data(freqsNLsgmplus)
#' 
#' # sample relatives of one profile
#' x1 <- sample.profiles(N=1,freqsNLsgmplus)
#' x1.sibs <- sample.relatives(x=x1,N=10^3,type="FS")
#' nrow(x1.sibs) # 10^3
#' 
#' # sample relatives of many profiles
#' x2 <- sample.profiles(N=10^3,freqsNLsgmplus)
#' x2.sibs <- sample.relatives(x=x2,N=1,type="FS")
#' nrow(x2.sibs) # 10^3
#' 
#' @export
sample.relatives <- function(x,N,type="FS",freqs=get.freqs(x),markers=names(freqs),theta=0){
    x.markers <- get.markers(x) # does a check on the column names of x as well
    freqs.markers <- names(freqs)
    
    # check if profiles and frequencies are available for these markers
    if (!all(markers %in% freqs.markers)){      
      stop("Allele frequencies unavailable for marker(s) ",paste(markers[!markers %in% freqs.markers],collapse=", "))}
    if (!all(markers %in% x.markers)){      
      stop("x does not contain marker(s) ",paste(markers[!markers %in% x.markers],collapse=", "))}
    
    ret <-  matrix(integer(),nrow=N*nrow(x),ncol=2*length(markers))
    colnames(ret) <- paste(rep(markers,each=2),1:2,sep=".")
    
    k <- ibdprobs(type)
    
    #first check whether x is one profile or many
    if (nrow(x)>1) {
      ismany <- TRUE
      # there are nrow(x) profiles, each gets N relatives
      # i.rel maps the nrow(x)*N relatives to the corresponding profile in x
      i.rel <- rep(1:nrow(x),each=N) # e.g. for N=2: 1,1,2,2,3,3,...
    }else{
      ismany <- FALSE
    } 
    
    for (m in markers){
      ind.ret <- match(m,markers)*2+c(-1,0)
      ind.x <- match(m,x.markers)*2+c(-1,0)    
      
      # look up allele freqs
      fr <- as.vector(freqs[[m]]);    fr.n <- length(fr)
      
      # decide which relatives have 0,1,2 alleles ibd with x
      ibd <- sample(0:2,size=N*nrow(x),replace=TRUE,prob=k) 
      which.0 <- which(ibd==0L); which.1 <- which(ibd==1L); which.2 <- which(ibd==2L)
      
      a <- as.integer(x[,ind.x[1]])
      b <- as.integer(x[,ind.x[2]])
      
      if (!ismany){
        # x is one profile, for which N relatives are sampled
        # three cases: the relatives have 0,1 or 2 alleles ibd
        
        if (!any(is.na(c(a,b)))){
          # 0 ibd: enumerate possible genotypes and compute their pr.'s
          G <- Zcomb.pairs(nn = length(fr))
          G.pr <- (2-(G[,1]==G[,2]))*Zprnextalleles(ij = G, seen = matrix(c(rep(a,nrow(G)),rep(b,nrow(G))),ncol=2),fr = fr,theta = theta)
          # sample from the joint dist of these genotypes
          G.ind <- sample.int(nrow(G),size=length(which.0),replace=TRUE,prob=G.pr)
          ret[which.0,ind.ret] <- G[G.ind,]
          
          # 1 ibd: one allele is taken from x, the other is sampled conditional on the alleles of x
          ret[which.1,ind.ret[1]] <- sample(c(a,b),length(which.1),replace=TRUE) #ibd allele
          # derive pr. distribution of the sampled allele
          A.pr <- pr.next.allele(seq_along(fr),seen = matrix(c(rep(a,fr.n),rep(b,fr.n)),ncol=2),fr = fr,theta = theta)
          ret[which.1,ind.ret[2]] <- sample.int(fr.n,size = length(which.1),replace = TRUE,prob = A.pr)  
          
          # 2 ibd: both allles from x
          ret[which.2,ind.ret[1]] <- a
          ret[which.2,ind.ret[2]] <- b
        }
      }else{
        # x contains multiple profiles; sample N relatives for each
        
        if (any(is.na(c(a,b)))) stop("nrow(x)>1 is not supported when x contains NAs")
        
        # 0 ibd: enumerate possible genotypes and compute their pr.'s
        G <- Zcomb.pairs(nn = length(fr))
        G.n <- nrow(G)
        # determine for each row of G the pr. dist of the next two alleles
        G.dist <- apply(G,1,function(ab) (2-(G[,1]==G[,2]))*Zprnextalleles(ij = G, seen = matrix(c(rep(ab[1],nrow(G)),rep(ab[2],nrow(G))),ncol=2),fr = fr,theta = theta))
        # determine for each row of x[irel[which.0],] which pr. distribution applies
        # match x to G, but remember that G is ordered st G[,1]>G[,2]
        swap <- b>a
        tmp <- b; b[swap] <- a[swap];a[swap] <- tmp[swap]
        
        a.0 <- a[i.rel[which.0]]; b.0 <- b[i.rel[which.0]]
        rel.G.i <- (fr.n*(b.0-1)-(b.0)*(b.0-1)/2)+a.0
        
        for(j in seq_len(G.n)){
          rel.G.ind <- which(rel.G.i==j)
          ret[which.0[rel.G.ind],ind.ret] <- G[sample.int(n = G.n,size = length(rel.G.ind),replace = TRUE,prob = G.dist[,j]),]        
        }
        
        # 1 ibd: one allele is ibd with x
        a.1 <- a[i.rel[which.1]]; b.1 <- b[i.rel[which.1]]      
        ret[which.1,ind.ret[1]] <- x[cbind(i.rel[which.1],sample(ind.x,size = length(which.1),replace = TRUE))]
        
        # sample the other allele
        # determine for each row of G the pr. dist of the next allele
        A.dist <- apply(G,1,function(ab) Zprnextalleles(ij = matrix(seq_len(fr.n),ncol=1), 
                                                                      seen = matrix(c(rep(ab[1],fr.n),rep(ab[2],fr.n)),ncol=2),
                                                                      fr = fr,theta = 0))  
        rel.A.i <- (fr.n*(b.1-1)-(b.1)*(b.1-1)/2)+a.1
        for(j in seq_len(G.n)){
          rel.A.ind <- which(rel.A.i==j)        
          ret[which.1[rel.A.ind],ind.ret[2]] <- sample.int(n = fr.n,size = length(rel.A.ind),replace = TRUE,prob = A.dist[,j])
        }
        
        # 2 ibd
        ret[which.2,ind.ret[1]] <- a[i.rel[which.2]]
        ret[which.2,ind.ret[2]] <- b[i.rel[which.2]]      
      }    
    }
    
    class(ret) <- c("profiles",class(ret))
    attr(ret,"freqs") <- freqs  
    ret
}
NULL
#' Sample random profile pairs with given relationship (sibs, parent/offspring, etc.)
#' 
#' @param N The number of pairs to be sampled (integer).
#' @param type A character string giving the type of relative. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param markers A character vector naming the markers of the resulting sample. Default to all markers of the \code{freqs} argument.
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @details The function randomly samples \eqn{N} pairs of DNA profiles according to the specified allelic frequencies. It returns two matrices containing profiles. The \eqn{i}'th profile in the first and the second matrix are sampled as relatives.
#' @return A list containing two integer matrices of class \code{profiles}:
#'          \enumerate{
#'            \item \code{x1} An integer matrix with \eqn{N} profiles.
#'            \item \code{x2} An integer matrix with \eqn{N} profiles.
#'          }
#'          
#'          
#' @seealso \code{\link{sample.profiles}}, \code{\link{sample.relatives}},\code{\link{ki}},\code{\link{ibs.pairs}}
#' @examples
#' ## Compare the number of IBS alleles of simulated parent/offspring pairs
#' ## with simulated unrelated pairs
#' 
#' data(freqsNLsgmplus)
#' 
#' #sample PO pairs and UN pairs
#' po.pairs <- sample.pairs(N=10^4,"PO",freqsNLsgmplus)
#' unr.pairs <- sample.pairs(N=10^4,"UN",freqsNLsgmplus)
#' 
#' #count the IBS alleles
#' po.pairs.ibs <- ibs.pairs(x1=po.pairs$x1,x2=po.pairs$x2)
#' unr.pairs.ibs <- ibs.pairs(x1=unr.pairs$x1,x2=unr.pairs$x2)
#' 
#' #plot together in a histogram
#' hist(po.pairs.ibs$ibs,breaks=0:20,xlim=c(0,20),
#' col="#FF0000FF",main="PO pairs vs. UN pairs",xlab="IBS")
#' hist(unr.pairs.ibs$ibs,breaks=0:20,col="#0000FFBB",add=TRUE)
#' @export
sample.pairs <- function(N=1,type="FS",freqs,markers=names(freqs)){
  k <- ibdprobs(type) # look up ibd probs
  
  #first sample N profiles, then sample the other N profiles of given type w.r.t. the first profiles
  prof1 <- sample.profiles(N=N,markers = markers,freqs = freqs)
  prof2 <- sample.relatives(x=prof1,N=1,markers = markers,type=k)
  list(x1=prof1,x2=prof2)
}
NULL