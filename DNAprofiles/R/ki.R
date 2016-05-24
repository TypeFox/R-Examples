#' Computes KIs for case profile(s) with all database profiles
#' 
#' @param x An integer matrix specifying either a single profile or a number of profiles. Alternatively an integer vector containing a single profile, e.g. obtained when a row is selected from a matrix of profiles.
#' @param db An integer matrix which is the database of profiles.
#' @param hyp.1 A character string giving the hypothesis in the numerator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param hyp.2 A character string giving the hypothesis in the denominator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param markers Character vector stating the markers to use in the KI computation. Defaults to the intersection of the markers of \code{x1} and \code{x2}.
#' @param theta Numeric value specifying the amount of background relatedness.
#' @param ret.per.marker Logical. If TRUE, return a matrix of KIs, where the columns correspond to markers.
#' @param precomputed.kis (optionally) a list of precomputed KIs, returned by \code{\link{ki.precompute}}. This speeds up the computation when multiple profiles are run against the db (i.e. \code{x} has more than one row).
#' @examples
#' 
#' data(freqsNLsgmplus)
#' fr <- freqsNLsgmplus
#'
#' # sample a profile, a database and compute the Sibling Index (SI) with all database members
#' x <- sample.profiles(N=1,fr)
#' db <- sample.profiles(N=10^4,fr)
#' si <- ki.db(x,db=db,"FS")
#'
#' # estimate the exceedance probabilities of an SI-threshold
#' t <- 1 # choose threshold SI=1
#' x <- sample.profiles(N=1,fr)
#' sibs <- sample.relatives(x,N=10^4,type="FS")
#' unrelated <- sample.profiles(N=10^4,fr)
#' mean(ki.db(x,db=sibs,"FS")>=t) # the vast majority of true siblings has an SI>=1
#' mean(ki.db(x,db=unrelated,"FS")>=t) # a few percent of unrelated persons have SI >= 1
#' 
#' # estimate distribution of SI for true siblings and unrelated persons
#' x <- sample.profiles(N=1,fr) #sample profile
#' sibs <- sample.relatives(x,N=10^4,type="FS") #sample sibs
#' unrelated <- sample.profiles(N=10^4,fr) #sample unrelated persons
#' 
#' sibs.si <- ki.db(x,db=sibs,"FS") #compute SI for true siblings
#' unrelated.si <- ki.db(x,db=unrelated,"FS") #compute SI for unrelated persons
#' #plot density estimates of SI
#' plot(density(log10(unrelated.si)),xlim=c(-10,10),lty="dashed",
#'     xlab=expression(log[10](SI)),main="SI for true sibs and unrelated profiles")
#' lines(density(log10(sibs.si)))
#' @export
ki.db <- function(x,db,hyp.1,hyp.2="UN",freqs=get.freqs(x),markers=intersect(get.markers(x),get.markers(db)),theta=0, ret.per.marker=FALSE,precomputed.kis){  
  x <- Zassure.matrix(x)
  db <- Zassure.matrix(db)
  
  if (missing(precomputed.kis)||(nrow(x)==1)){
    if (!identical(ibdprobs(hyp.2),ibdprobs("UN"))){
      ki.db(x = x,db = db,hyp.1 = hyp.1,hyp.2 = "UN",freqs = freqs, markers = markers, theta = theta, ret.per.marker = ret.per.marker)/
        ki.db(x = x,db = db,hyp.1 = hyp.2,hyp.2 = "UN",freqs = freqs, markers = markers, theta = theta, ret.per.marker = ret.per.marker)        
    }    
    
    # check if freqs available for markers
    freqs.markers <- names(freqs)
    if (!all(markers %in% freqs.markers)){      
      stop("Allele frequencies unavailable for marker(s) ",paste(markers[!markers %in% freqs.markers],collapse=", "))}
    
    # check for off-ladder alleles that would lead to a crash
    min.x <- min(x,na.rm = TRUE)
    max.x <- max(x,na.rm = TRUE)
    if (min.x<1L) stop("Alleles should be positive integers")
    freqs.L <- sapply(freqs[markers],length)
    if (max.x>max(freqs.L)) stop("x contains allele that is not in freqs")
    k <- ibdprobs(hyp.1)
    fr.mat <- Zfreqs.to.mat(freqs = freqs[markers])
  }
  
  # check if all markers of x1 and x2 are present
  x.markers <- get.markers(x)
  db.markers <- get.markers(db)
  
  # check if profiles are available for these markers
  if (!all(markers %in% x.markers)){      
    stop("x1 does not contain marker(s) ",paste(markers[!markers %in% x.markers],collapse=", "))}
  if (!all(markers %in% db.markers)){      
    stop("x2 does not contain marker(s) ",paste(markers[!markers %in% db.markers],collapse=", "))}
 
  
  # actual ki computation below
  x.ind <- match(markers,x.markers)-1
  db.ind <- match(markers,db.markers)-1
  
  if (nrow(x)==1){
    ## single profile vs db
    ret <- Zki(x1 = x,manytomany = FALSE,x2 = db,x1ind = x.ind,x2ind = db.ind,fr = fr.mat,
                             k0 = k[1], k1 = k[2],k2 = k[3],theta = theta,retpermarker = as.integer(ret.per.marker))
  }else{
    ## we run multiple (k) profiles against the database (N)
    # we return a N*k matrix with KIs
    
    if (ret.per.marker) stop("Can not return per marker when multiple profiles are searched against db")
    
    # if the user supplies a lookup table of KIs, then use this,
    # else just compute the KIs
    if (!missing(precomputed.kis)){
      nm <- paste(rep(markers,each=2),1:2,sep=".")
      if (any(!markers%in%names(precomputed.kis))) stop("Precomputed KIs not available for marker(s)",
                                                        paste(markers[!markers%in%names(precomputed.kis)],sep=", "))
      ret <- ZcompKItargetsdbwithtable(precomputed.kis[markers],x[,nm,drop=FALSE],db[,nm,drop=FALSE])
    }else{
      ret <- matrix(NA,nrow = nrow(db),ncol = nrow(x))
      for(i in seq_len(nrow(x))){
        ret[,i] <- Zki(x1 = x[i,,drop=FALSE],manytomany = FALSE,x2 = db,x1ind = x.ind,x2ind = db.ind,fr = fr.mat,
                                     k0 = k[1], k1 = k[2],k2 = k[3],theta = theta,retpermarker = as.integer(ret.per.marker)) }
    }
  }
  
  if (ret.per.marker) colnames(ret) <- markers
  ret
}
NULL
#' Computes Kinship Indices (KIs) for pairs of profiles
#' 
#' @param x1 An integer matrix with \eqn{N} profiles.
#' @param x2 An integer matrix with \eqn{N} profiles.
#' @param hyp.1 A character vector giving the hypothesis in the numerator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param hyp.2 A character vector giving the hypothesis in the denominator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each marker, named after that marker. 
#' @param markers Character vector stating the markers to use in the KI computation. Defaults to the intersection of the markers of \code{x1} and \code{x2}.
#' @param theta Numeric value specifying the amount of background relatedness.
#' @param ret.per.marker Logical. If TRUE, return a matrix of KIs, where the columns correspond to markers.
#' @seealso \link{ibs.pairs}
#' @export
#' @examples
#' data(freqsNLngm)
#' fr <- freqsNLngm
#' sibs1 <- sample.profiles(1e3,fr) # sample profiles
#' sibs2 <- sample.relatives(sibs1,1,type="FS",freqs=fr) #sample 1 sib for each profile
#' #compute ki for all pairs
#' ki(sibs1,sibs2,hyp.1="FS",hyp.2="UN")
#' @export
ki <- function(x1,x2,hyp.1,hyp.2="UN",freqs=get.freqs(x1),markers=intersect(get.markers(x1),get.markers(x2)),theta=0, ret.per.marker=FALSE){
  if (!identical(ibdprobs(hyp.2),ibdprobs("UN"))){
    return(ki(x1 = x1,x2 = x2,hyp.1 = hyp.1,hyp.2 = "UN",freqs = freqs,markers = markers,theta = theta,ret.per.marker = ret.per.marker)/
           ki(x1 = x1,x2 = x2,hyp.1 = hyp.2,hyp.2 = "UN",freqs = freqs,markers = markers,theta = theta,ret.per.marker = ret.per.marker))
  }
  
  x1 <- Zassure.matrix(x1)
  x2 <- Zassure.matrix(x2)
  if (nrow(x1)!=nrow(x2)) stop("Number of profiles in x1 is unequal to number of profiles in x2")
  
  # first check if all markers of x1 and x2 are present and allele ladders are available  
  x1.markers <- get.markers(x1)
  x2.markers <- get.markers(x2)
  freqs.markers <- names(freqs)
  # check if profiles and frequencies are available for these markers
  if (!all(markers %in% x1.markers)){      
    stop("x1 does not contain marker(s) ",paste(markers[!markers %in% x1.markers],collapse=", "))}
  if (!all(markers %in% x2.markers)){      
    stop("x2 does not contain marker(s) ",paste(markers[!markers %in% x2.markers],collapse=", "))}
  if (!all(markers %in% freqs.markers)){      
    stop("Allele frequencies unavailable for marker(s) ",paste(markers[!markers %in% freqs.markers],collapse=", "))}
  
  # check for off-ladder alleles that would lead to a crash
  min.x1 <- min(x1,na.rm = TRUE)
  max.x1 <- max(x2,na.rm = TRUE)
  if (min.x1<1L) stop("Alleles should be positive integers")
  freqs.L <- sapply(freqs[markers],length)
  if (max.x1>max(freqs.L)) stop("x1 contains allele that is not in freqs")
  
  # actual ki computation
  k <- ibdprobs(hyp.1)
  x1.ind <- match(markers,x1.markers)-1
  x2.ind <- match(markers,x2.markers)-1
  fr.mat <- Zfreqs.to.mat(freqs = freqs[markers])
  ret <- Zki(x1 = x1,manytomany = TRUE,x2 = x2,x1ind = x1.ind,x2ind = x2.ind,fr = fr.mat,
                    k0 = k[1], k1 = k[2],k2 = k[3],theta = theta,retpermarker = as.integer(ret.per.marker))
  
  if (ret.per.marker) colnames(ret) <- markers
  ret
}
NULL