ZfindInterval <- function(x, vec, rightmost.closed = FALSE, all.inside = FALSE)
{
  #if(any(is.na(vec)))
  #  stop("'vec' contains NAs")
  #if(is.unsorted(vec))
  #  stop("'vec' must be sorted non-decreasingly")
  stopifnot(!all.inside)
  ZfindIntervalcpp(x,vec)
  #.Internal(findInterval(as.double(vec), as.double(x),rightmost.closed, all.inside))
}

Zcutright.str <- function(x,n){
  #cuts off n+1 characters from the right
  substr(x,start=1,stop=(nchar(x)-n))
}

Zprofile.names <- function(profile){
  # extracts the locus names from a profile
  # these are either colnames when profile is a 1-row matrix (e.g. as a result from sample.profiles(N=1,..))
  # or names when profile is a vector (e.g. when db[1,,drop=TRUE] is used)
  if (is.null(colnames(profile))) ret <- names(profile)
  else{
    ret <- colnames(profile)
  } 
  ret
}

Znames.to.loci <- function(x){
  # convert colnames of db/target to character vectors of loci
  # e.g. D2S1338.1 D2S1338.2 D3S1358.1 D3S1358.2     FGA.1     FGA.2
  # to c("D2S1338", "D3S1358",...)
  if (length(x)==0) stop("Can not read loci from column names of profile(s).")
  as.vector(sapply(x[seq(from=1,to=length(x),by=2)],function(y) Zcutright.str(y,2)) )
}

Zloci <- function(x){
  # extract the names of the loci from a profiles object
  Znames.to.loci(Zprofile.names(x))
}

Zassure.matrix <- function(x){
  # sometimes x is a db (matrix)
  # sometimes it is a single profile (vector) e.g. when x <- db[1,] without drop=FALSE
  # this functions takes x and makes it an n x (2*nloci) matrix (possibly n=1)
  if (!is.matrix(x)){
    tmp <- attr(x,"freqs")
    x <- matrix(x,nrow=1,dimnames=list("",Zprofile.names(x)))
    attr(x,"freqs") <- tmp
  }
  x
}

Zall.unordered.pairs.int <- function(nn){ # order is wrong, will be removed
  # all genotypes at a locus, second index varies fastest
  cbind(rep(1:nn,times=(nn:1)),unlist(sapply(1:nn,function(n) n:nn)))
}

Zcomb.pairs <- function(nn,firstindexfastest=TRUE){
  # all genotypes at a locus, first index varies fastest
  if (firstindexfastest){
    return(cbind(unlist(sapply(1:nn,function(n) n:nn)),rep(1:nn,times=(nn:1))))
  }else{
    return(cbind(rep(1:nn,times=(nn:1)),unlist(sapply(1:nn,function(n) n:nn))))
  }
}

#Zchecktype <- function(type){
#  # check if the type of relative is known, i.e. is one of FS, PO, PO, UN, HS, AV, FC, SC
#  sapply(type,function(type) if (is.null(ibdprobs[[type]])) stop("Unknown type of relative: ",type,". Choose one of ",paste(names(ibdprobs),collapse=", ")))
#}

Zibdpr <- list() # prob. of sharing 0,1,2 alleles ibd
Zibdpr[["ID"]] <- c(0,0,1)
Zibdpr[["FS"]] <- c(1/4,1/2,1/4)
Zibdpr[["PO"]] <- Zibdpr[["PC"]] <- c(0,1,0)
Zibdpr[["UN"]] <- c(1,0,0)
Zibdpr[["HS"]] <- c(1/2,1/2,0)
Zibdpr[["AV"]] <- c(1/2,1/2,0)
Zibdpr[["FC"]] <- c(3/4,1/4,0)
Zibdpr[["SC"]] <- c(15/16,1/16,0)

Zchecktheta <- function(theta) if (!((theta>=0)&(theta<=1))) stop("theta must be non-negative and at most 1!")

Zdiststomatrix.X <- function(dists){
  X.n <- sapply(dists,function(x) length(x$x))
  X <- matrix(NA,nrow=max(X.n),ncol=length(dists))  
  for(j in seq_along(dists))    X[1:X.n[j],j] <- dists[[j]]$x
  X
}
Zdiststomatrix.P <- function(dists){
  P.n <- sapply(dists,function(x) length(x$fx))
  P <- matrix(NA,nrow=max(P.n),ncol=length(dists))  
  for(j in seq_along(dists))    P[1:P.n[j],j] <- dists[[j]]$fx
  P
}

Zfind.subsets.with.max.product <- function(x,max.product){
  # function selects subsets of x s.t. the product of their elements is smaller than max.product
  # e.g. x=c(1,2,3,4), max.product=10 yields list(c(1,2,3), c(4))
  ret <- list()
  for (i in seq_along(x)){
    if (i==1) ret[[i]] <- i
    else{      
      if (prod(c(x[ret[[length(ret)]]],x[i]))<max.product){
        ret[[length(ret)]] <- c(ret[[length(ret)]],i)
      }else{
        ret[[length(ret)+1]] <- i
      }
    }
  }
  ret
}
NULL
#' Normalizes allele frequencies such that their sum is 1
#' 
#' @param freqs list of per locus allele frequencies
#' @return list
#' @examples 
#' 
#' data(freqsNLsgmplus)
#' fr0 <- normalize.freqs(freqsNLsgmplus)
#' stopifnot(all.equal(sapply(fr0,sum),setNames(rep(1,length(fr0)),names(fr0))))
#' @export
normalize.freqs <- function(freqs) lapply(freqs,function(x) x/sum(x))