#' Computes distribution of KI for profiles with stated relationship
#' 
#' Computes, per marker, the distribution of a Kinship Index (KI) comparing hypotheses \code{hyp.1} vs \code{hyp.2} for profiles with a given relationship (\code{hyp.true}). Optionally with respect to a profile \code{x}.
#' 
#' @param x (optional) An integer matrix specifying a single profile.
#' @param hyp.1 A character string giving the hypothesis in the numerator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param hyp.2 A character string giving the hypothesis in the denominator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param hyp.true A character string specifying the true relationship between the case profile and the other profile. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs.ki A list specifying the allelic frequencies that are used when computing the \eqn{KI}.
#' @param freqs.true (optionally) A list specifying the allelic frequencies that are used for computing the probabily distribution of the \eqn{KI} under \code{hyp.true}. When not provided, the function will use \code{freqs}. One might use different allelic frequencies \code{freqs.rel} when for example the case profile and relative come from some population, while \eqn{KI}s are computed with frequencies from another population.
#' @param markers Character vector stating the markers for which the KI distribution is derived. Defaults to the intersection of the markers of \code{freqs.ki} and \code{freqs.true}.
#' @param theta.ki numeric value specifying the amount of background relatedness.
#' @param theta.true numeric value specifying the amount of background relatedness.
#' @param min.freq Alleles with a frequency in \code{freqs.ki} and \code{freqs.true} smaller than this value will be set to frequency 0 to avoid numerical problems. Defaults to \code{.Machine$double.eps}, which is normally \code{2.220446e-16}.
#' @return A list of distributions, where a distribution is specified by a list with vectors \code{x}, \code{fx}.
#' @export
ki.dist <- function(x,hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki=get.freqs(x),freqs.true=freqs.ki,markers=intersect(names(freqs.ki),names(freqs.true)),theta.ki=0,theta.true=theta.ki,min.freq = .Machine$double.eps){ 
  if (missing(hyp.1)) stop("hyp.1 is missing")
  if (missing(x)){
    # unconditional (= for two profiles) ki dist
    ret <- sapply(markers,    function(m){
      freqs.ki[[m]][freqs.ki[[m]]<min.freq] <- 0
      freqs.true[[m]][freqs.true[[m]]<min.freq] <- 0
      y <- Zki.ibs.joint.dist.marker(k.hyp.1 = ibdprobs(hyp.1),k.hyp.2 = ibdprobs(hyp.2),k.hyp.true = ibdprobs(hyp.true),fr.ki = freqs.ki[[m]],fr.true = freqs.true[[m]],theta.ki = theta.ki,theta.true = theta.true)
      dist.unique.events(list(x=y$ki,fx= y$fx))},
      simplify = FALSE, USE.NAMES=TRUE)    
  }else{
    # ki dist of profile x and some profiles y, related to x by hyp.true
    dist <- ki.ibs.joint.dist(x=x,hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,freqs.ki=freqs.ki,freqs.true=freqs.true,markers = markers,theta.ki=theta.ki,theta.true=theta.true)
    ret <- lapply(dist, function(y) dist.unique.events(list(x=y$ki,fx=y$fx)))    
  }
  ret
}
NULL
#' Computes distribution of number of IBS alleles for profiles with stated relationship
#' 
#' Computes, per locus, the distribution of a Kinship Index (KI) comparing hypotheses \code{hyp.1} vs \code{hyp.2} for profiles with a given relationship (\code{hyp.true}), optionally with respect to the case profile (e.g. \code{"FS"} for full siblings).
#' 
#' @param x (optional) An integer matrix specifying a single profile.
#' @param hyp.true A character string specifying the true relationship between the two profiles. Forwarded to \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs A list specifying the allelic frequencies.
#' @param markers A character vector giving the markers for which the distribution is derived. Default to the markers of \code{freqs}.
#' @param theta numeric value specifying the amount of background relatedness.
#' @return A list of distributions, where a distribution is specified by a list with vectors \code{x}, \code{fx}.
#' @export
ibs.dist <- function(x,hyp.true="UN",freqs=get.freqs(x),markers=names(freqs),theta=0){
  #  ibs dist at all loci
  jd <- ki.ibs.joint.dist(x,hyp.1="UN",hyp.2="UN",hyp.true=hyp.true,markers=markers,freqs.ki=freqs,theta.true=theta)
  lapply(jd,function(y) dist.unique.events(list(x=y$ibs,fx=y$fx)) )
}
NULL
#' Computes joint distribution of KI and IBS for profiles with stated relationship
#' 
#' Computes, per locus, the joint distribution of the number of IBS alleles and a Kinship Index (KI) comparing hypotheses \code{hyp.1} vs \code{hyp.2} for profiles with a given relationship (\code{hyp.true}), optionally with respect to the case profile (e.g. \code{"FS"} for full siblings).
#' 
#' @param x (optional) An integer matrix specifying a single profile.
#' @param hyp.1 A character string giving the hypothesis in the numerator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param hyp.2 A character string giving the hypothesis in the denominator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param hyp.true A character string specifying the true relationship between the case profile and the other profile. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs.ki A list specifying the allelic frequencies that are used when computing the \eqn{KI}.
#' @param freqs.true (optionally) A list specifying the allelic frequencies that are used for computing the probabily distribution of the \eqn{KI} under \code{hyp.true}. When not provided, the function will use \code{freqs}. One might use different allelic frequencies \code{freqs.rel} when for example the case profile and relative come from some population, while \eqn{KI}s are computed with frequencies from another population.
#' @param markers Character vector stating the markers for which the KI distribution is derived. Default to the intersection of the markers of \code{freqs.ki} and \code{freqs.true}.
#' @param theta.ki numeric value specifying the amount of background relatedness.
#' @param theta.true numeric value specifying the amount of background relatedness.
#' @return A list of distributions, where a distribution is specified by a list with vectors \code{ki}, \code{ibs}, \code{fx}.
#' @export
ki.ibs.joint.dist <- function(x,hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki=get.freqs(x),freqs.true=freqs.ki,markers=intersect(names(freqs.ki),names(freqs.true)),theta.ki=0,theta.true=theta.ki){
  if (!all(markers %in% names(freqs.ki))){
    stop("Freqs.ki does not contain marker(s) ",paste(markers[!markers %in% names(freqs.ki)],collapse=", "))}
  if (!all(markers %in% names(freqs.true))){
    stop("Freqs.true does not contain marker(s) ",paste(markers[!markers %in% names(freqs.ki)],collapse=", "))}
  
  if (missing(x)){
    # unconditional ki,ibs joint dist
    ret <- sapply(markers, function(L) Zki.ibs.joint.dist.marker(k.hyp.1 = ibdprobs(hyp.1),k.hyp.2 = ibdprobs(hyp.2),
                                                                 k.hyp.true = ibdprobs(hyp.true),
                                                          fr.ki = freqs.ki[[L]],fr.true = freqs.true[[L]],
                                                          theta.ki = theta.ki,theta.true = theta.true),
                  simplify=FALSE,USE.NAMES=TRUE)
  }else{
    
    x <- Zassure.matrix(x)
    x.markers <- get.markers(x)
    if (!all(markers%in%x.markers)){
      stop("x does not contain marker(s) ",paste(markers[!markers %in% x.markers],collapse=", "))}
    
    if (nrow(x)>1) warning("nrow(x)>1, only first profile is used!")
    
    Zchecktheta(theta.ki);Zchecktheta(theta.true)  
    
    ret <- list()
    for(m in markers){
      tmp <- Zcond.ki.ibs.joint.dist.marker(a = as.integer(x[,paste(m,".1",sep="")]),b = as.integer(x[,paste(m,".2",sep="")]),
                                            k.hyp.1 = ibdprobs(hyp.1),k.hyp.2 = ibdprobs(hyp.2),k.hyp.true = ibdprobs(hyp.true),
                                            fr.ki = freqs.ki[[m]],fr.true = freqs.true[[m]],theta.ki = theta.ki,theta.true = theta.true)
      ret[[m]] <- list(fx=tmp[,1], ki = tmp[,2], ibs = tmp[,3])  
    }
    
  }
  ret
}
NULL
Zcond.ki.ibs.joint.dist.marker <-function(a,b,k.hyp.1,k.hyp.2=c(1,0,0),k.hyp.true=c(1,0,0),fr.ki,fr.true=fr.ki,theta.ki=0,theta.true=theta.ki){
  # function derives the joint dist of ki,ibs at a marker for a fixed profile with alleles a and b
    
  if (!any(is.na(c(a,b)))){
    # two cases: homozygous or heterozygous
    if (a==b){
      # homozygous: enumerate all combinations: (aa,zz), (aa,az), (aa,aa) (consider 2 alleles)
      fr0.ki <- c(fr.ki[a],  1-sum(fr.ki[a]))
      fr0.true <- c(fr.true[a],  1-sum(fr.true[a]))
      # enumerate profiles
      x1 <- matrix(c(1L,1L),nrow=1)
      x2 <- Zcomb.pairs(2L,firstindexfastest=FALSE)
      
            # compute possible KIs
      lr <- Zki(x1 = x1,manytomany = FALSE,x2 = x2,x1ind = 0,x2ind = 0,fr = matrix(fr0.ki,ncol=1),
                              k0 = k.hyp.1[1],k1 = k.hyp.1[2],k2 = k.hyp.1[3],theta = theta.ki,retpermarker = FALSE)/
        Zki(x1 = x1,manytomany = FALSE,x2 = x2,x1ind = 0,x2ind = 0,fr = matrix(fr0.ki,ncol=1),
                          k0 = k.hyp.2[1],k1 = k.hyp.2[2],k2 = k.hyp.2[3],theta = theta.ki,retpermarker = FALSE)
      
      ibs <- c(2L,1L,0L)  
      
      # and the probability by which these occur. use that Pr(LR=z|Htrue) = LR_{Htrue,Hunr} Pr(LR=z|Hunr)
      lr.hyp.true <- Zki(x1 = x1,manytomany = FALSE,x2 = x2,x1ind = 0,x2ind = 0,fr = matrix(fr0.true,ncol=1),
                                       k0 = k.hyp.true[1],k1 = k.hyp.true[2],k2 = k.hyp.true[3],theta = theta.true,retpermarker = FALSE)
      p.hyp.unr <- (2-(x2[,1]==x2[,2]))*Zprnextalleles(ij = x2,seen = matrix(x1,nrow=nrow(x2),ncol=ncol(x2),byrow=TRUE),fr = fr0.true,theta = theta.true)
      
      p.x <- lr.hyp.true*p.hyp.unr
    }else{
      # heterozygous. enumerate all combinations: (ab,zz), (ab,az), ..., (ab,ab) (consider 3 alleles)
      fr0.ki <- c(fr.ki[c(a,b)],  1-sum(fr.ki[c(a,b)]))
      fr0.true <- c(fr.true[c(a,b)],  1-sum(fr.true[c(a,b)]))
      # enumerate profiles
      x1 <- matrix(1:2,nrow=1)
      x2 <- Zcomb.pairs(3L,firstindexfastest=FALSE)
      
      # compute possible KIs
      lr <- Zki(x1 = x1,manytomany = FALSE,x2 = x2,x1ind = 0,x2ind = 0,fr = matrix(fr0.ki,ncol=1),
                              k0 = k.hyp.1[1],k1 = k.hyp.1[2],k2 = k.hyp.1[3],theta = theta.ki,retpermarker = FALSE)/
            Zki(x1 = x1,manytomany = FALSE,x2 = x2,x1ind = 0,x2ind = 0,fr = matrix(fr0.ki,ncol=1),
                          k0 = k.hyp.2[1],k1 = k.hyp.2[2],k2 = k.hyp.2[3],theta = theta.ki,retpermarker = FALSE)

      ibs <- c(1L,2L,1L,1L,1L,0L)  
      
      # and the probability by which these occur. use that Pr(LR=z|Htrue) = LR_{Htrue,Hunr} Pr(LR=z|Hunr)
      lr.hyp.true <- Zki(x1 = x1,manytomany = FALSE,x2 = x2,x1ind = 0,x2ind = 0,fr = matrix(fr0.true,ncol=1),
                                       k0 = k.hyp.true[1],k1 = k.hyp.true[2],k2 = k.hyp.true[3],theta = theta.true,retpermarker = FALSE)
      p.hyp.unr <- (2-(x2[,1]==x2[,2]))*Zprnextalleles(ij = x2,seen =  matrix(x1,nrow=nrow(x2),ncol=ncol(x2),byrow=TRUE),
                                                       fr = fr0.true,theta = theta.true)
      
      p.x <- lr.hyp.true*p.hyp.unr  
    }
  }else{
    # simply return 1
    p.x=1 ; ibs=0; lr=1;
  }
  
  p.notna <- (!is.na(p.x))
  p.nonzero <- p.x[p.notna]>0 # only retain the events with non-zero probability
  
  p.x <- p.x[p.notna][p.nonzero]
  lr <- lr[p.notna][p.nonzero]
  ibs <- ibs[p.notna][p.nonzero]
  
  if (any(is.na(lr))) stop("NA or NaNs encountered")
  cbind(p.x,lr,ibs)
}

Zki.ibs.joint.dist.marker <- function(k.hyp.1,k.hyp.2=c(1,0,0),k.hyp.true=c(1,0,0),fr.ki,fr.true=fr.ki,theta.ki=0,theta.true=0){
  # determines the unconditional ki,ibs joint dist for one marker
  
  # first determine all genotypes with fr.
  G <- enum.profiles(list(locus=fr.true))
  G.pr <- rmp(G,theta = theta.true)
  
  # then determine ki dist for all genotypes
  X <- list()
  for(i in seq_along(G[,1])){
    if (G.pr[i]>0) X[[i]] <- cbind(Zcond.ki.ibs.joint.dist.marker(as.integer(G[i,1]),as.integer(G[i,2]),k.hyp.1=k.hyp.1,k.hyp.2=k.hyp.2,k.hyp.true=k.hyp.true,fr.ki=fr.ki,fr.true=fr.true,theta.ki=theta.ki,theta.true=theta.true),G.pr[i])
  }
  X.matrix <- do.call(rbind,X) 
  # the columns are respectively pr(ki|x), ki, ibs, pr(x)
  # we are interested in the ki,ibs joint dist, so we multiply column 1 and 4
  list(fx=X.matrix[,1]*X.matrix[,4],ki=X.matrix[,2],ibs=X.matrix[,3]) # events are not unique!
}