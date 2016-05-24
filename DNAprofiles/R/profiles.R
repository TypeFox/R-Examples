#' @name print.profiles
#' @title print profiles object
#' @method print profiles
#' @param x profiles object
#' @param ... passed on to \code{print}
#' @aliases print.profiles
#' @description Profiles are stored in a profiles object, which is merely an integer matrix together with allelic frequencies stored as an attribute "freqs".
#' @examples data(freqsNLsgmplus)
#'            x<- sample.profiles(1,freqsNLsgmplus)
#'            print(x)
#' @export
print.profiles <- function(x,...){
   tmp <- attr(x,"freqs")
   attr(x,"freqs") <- NULL
   print(unclass(x))
   attr(x,"freqs") <- tmp
   class(x) <- c("profiles",class(x))
   invisible(x)
 }
NULL
#' @name profiles
#' @title profiles object
#' @param x list of character vectos with genotypes
#' @param labels allele levels per locus
#' @description Profiles are stored in a profiles object, which is merely an integer matrix. The profiles function creates such a matrix from the genotypes of profiles stored in a list. Alleles are encoded as integers following the order in the labels argument. See the example for details.
#' @examples data(freqsNLngm)
#'x <- profiles(list(D1S1656="12/15",D2S441="11/14",D2S1338="19/19",D3S1358="16/18",
#'                    FGA="22/24",D8S1179="11/12",D10S1248="13/14",TH01="7/9",VWA="16/17",
#'                    D12S391="17/20",D16S539="12/13",D18S51="13/17",
#'                    D19S433="13/14",D21S11="27/31.2",D22S1045="15/16"),
#'              labels = get.labels(freqs = freqsNLngm))
#' x
#' @export
profiles <- function(x,labels){
  nm <- names(x)
  
  x.len <- sapply(x,length)
  n <- max(x.len)
  
  if (all(substr(nm,start = nchar(nm)-2,stop = nchar(nm))%in%c(".1",".2"))){
    stop("handling of genotypes stored as one list item per marker is not implemented yet")
  }else{
    markers <- nm  
    if (!all(x.len==n)) stop("All elements of x should have equal length")
    
    ret <- matrix(integer(),nrow = n,ncol = 2*length(markers))
    colnames(ret) <- paste(rep(markers,each=2),1:2,sep=".")
    for(m in markers){
      x0 <- do.call(rbind,strsplit(x[[m]],"/"))
      ret[,paste(m,"1",sep=".")] <- match(x0[,1],table = labels[[m]])
      ret[,paste(m,"2",sep=".")] <- match(x0[,2],table = labels[[m]])    
    }  
  }
  ret
}
NULL
#' @name get.freqs
#' @title Retrieve the allele frequencies of a profiles object
#' @param x profiles object
#' @description A \code{\link{profiles}} object is an integer matrix with two columns per marker. This function extracts the allele frequencies from the attributes of the profiles object.
#' @examples data(freqsNLsgmplus)
#'            x<- sample.profiles(1,freqsNLsgmplus)
#'            stopifnot(identical(get.freqs(x),freqsNLsgmplus))
#' @export
get.freqs <- function(x){
  if (is.null(attr(x,"freqs"))) stop("x does not have freqs attribute")
  attr(x,"freqs")
}
#' @name get.labels
#' @title Retrieve the allele labels from the allele frequencies
#' @param freqs list of allele frequencies
#' @description A \code{\link{profiles}} object is an integer matrix with two columns per marker and may carry allele frequencies. This function extracts the labels from the allele frequencies.
#' @examples data(freqsNLsgmplus)
#'            x<- sample.profiles(1,freqsNLsgmplus)
#'            labels <- get.labels(freqsNLsgmplus)
#'            
#'            stopifnot(identical(labels,get.labels(get.freqs(x))))
#' @export
get.labels <- function(freqs) sapply(freqs,names,simplify = FALSE)
NULL
#' @name get.markers
#' @title Retrieve the markers of a profiles object
#' @param x profiles object
#' @description A \code{\link{profiles}} object is an integer matrix with two columns per marker. This function extracts the marker names from the column names of the profiles object.
#' @examples data(freqsNLsgmplus)
#'            x<- sample.profiles(1,freqsNLsgmplus)
#'            stopifnot(identical(get.markers(x),names(freqsNLsgmplus)))
#' @export
get.markers <- function(x){
  # obtain colnames (or names if x became a vector by subsetting)
  nm <- Zprofile.names(x)
  
  if (is.null(nm)) stop("Profiles object does not have column names")
  markers <- Zcutright.str(nm,2)
  # do we have every marker twice (consecutive)?
  
  if (!identical(markers[seq(from=1L,to=length(markers),by=2L)],
                 markers[seq(from=2L,to=length(markers),by=2L)])){
    stop("Each two consecutive column names of a profiles object should correspond to one marker")
  }
  
  if (!all(substr(nm,start = nchar(nm)-1,stop = nchar(nm))==c(".1",".2"))){
    stop("Each two consecutive column names of a profiles object should correspond to one marker")
  }
  
  markers[seq(from=1L,to=length(markers),by=2L)]
}
NULL
#' Obtain STR repeat numbers of profiles as character matrix
#'
#' @param x profiles object
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param two.cols.per.locus when FALSE (default), returns "a/b", else "a" and "b" in separate columns
#' @details Profiles are stored as an integer matrix, with the integers corresponding to repeat numbers found in the names attribute of the list of allelic frequencies. The current function converts the integer matrix to a character matrix with alleles.
#' @return A character matrix with a column for each locus.
#' @examples
#' data(freqsNLsgmplus)
#' profiles.to.chars(sample.profiles(N=2,freqs=freqsNLsgmplus))
#' @export
profiles.to.chars <- function(x,freqs=get.freqs(x),two.cols.per.locus=FALSE){
  LL <- Zcutright.str(colnames(x)[seq(to=length(colnames(x)),by=2)],2) #loci
  
  if (two.cols.per.locus){
    ret <- matrix(character(),nrow=nrow(x),ncol=(ncol(x)))
    for(i in seq_len(ncol(x)/2)){
      L <- LL[i]
      ret[,2*i-1] <- names(freqs[[L]])[x[,2*i-1]]
      ret[,2*i] <- names(freqs[[L]])[x[,2*i]]      
    }
    colnames(ret) <- colnames(x)
  }else{
    ret <- matrix(character(),nrow=nrow(x),ncol=(ncol(x)/2))
    for(i in seq_along(LL)){
      L <- LL[i]
      ret[,i] <- paste(names(freqs[[L]])[x[,2*i-1]],  names(freqs[[L]])[x[,2*i]],sep="/")    
    }
    colnames(ret) <- LL    
  }
  ret
}
NULL
#' Enumerate all attainable genotypes
#'
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @details Profiles are stored as an integer matrix, with the integers corresponding to repeat numbers found in the names attribute of the list of allelic frequencies. The current function converts the integer matrix to a character matrix with alleles.
#' @return A profiles object.
#' @examples
#' data(freqsNLsgmplus)
#' enum.profiles(freqsNLsgmplus[1:2])
#' @export
enum.profiles <- function(freqs){
  ret <- matrix(integer())
  
  fr.todo <- freqs[rev(seq_along(freqs))]
  while (length(fr.todo)>0){
    f <- fr.todo[[1]]
    A <- length(f) # number of alleles at locus
    g <- Zcomb.pairs(A,firstindexfastest=FALSE)
    colnames(g) <- paste(names(fr.todo)[1],c(1,2),sep=".")
    
    # second row index varies fastest
    i2 <- rep(seq_len(nrow(ret)),nrow(g))
    i1 <- rep(seq_len(nrow(g)),each=max(1,nrow(ret)))  
    
    ret <- cbind(g[i1,],  ret[i2,])    
    fr.todo <- fr.todo[-1]
  }
  
  class(ret) <- c("profiles",class(ret))
  attr(ret,which="freqs") <- freqs
  ret
}