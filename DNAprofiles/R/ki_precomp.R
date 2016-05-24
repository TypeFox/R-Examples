#' Pre-computes KIs for use with \code{ki.db} function
#' 
#' @param type A character string giving the type of KI. See \link{ibdprobs}.
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param theta numeric value specifying the amount of background relatedness.
#' @param markers Character vector stating the markers to use in the KI computation. Defaults to all markers contained in \code{freqs}.
#' @details In large scale simulation studies, it is sometimes useful to precompute KIs to speedup computations.
#' @return list A list of numeric vectors containing the KIs for all genotypic combinations at each locus.
#' @seealso \link{ki.db}
#' @export
#' @examples
#' \dontrun{
#' data(freqsNLngm);
#' n <- 1e6
#'
#' targets <- sample.profiles(N = 1e2,freqs = freqsNLngm)
#' db <- sample.profiles(N = 1e5,freqs = freqsNLngm)
#'
#' precomp <- ki.precompute(type = "FS",freqs = freqsNLngm)
#'
#' R1 <- ki.db(x = targets,db = db,hyp.1 = "FS")
#' R2 <- ki.db(x = targets,db = db,precomputed.kis = precomp) # a little faster
#'
#' all.equal(R1,R2)
#'
#'}
#' @export
ki.precompute <- function(type,freqs,markers=names(freqs),theta=0){
  sapply(markers,function(marker) Zprecompute.lrs.locus(marker,ibdprobs(type),freqs,theta=theta),
         simplify=FALSE,USE.NAMES=TRUE)
}
NULL
# the following function computes the KIs for all database geno's, conditional on a profile x
Zprecompute.lrs.locus.for.x <- function(x,locus,ki.type,fr,theta=0){
  L <- length(fr[[locus]])
  # all possible geno's, also reversed.. might optimize this at some point
  G <- cbind(rep(1:L,L),rep(1:L,each=L))
  colnames(G) <- paste(locus,c(".1",".2"),sep="")
  matrix(ki.db(x[,colnames(G)],db=G,hyp.1=ki.type,freqs=fr,theta=theta),nrow=L)
}

Zprecompute.lrs.for.x <- function(x,ki.type,fr,theta=0){
  sapply(Zloci(x),function(locus) Zprecompute.lrs.locus.for.x(x,locus,ki.type,fr,theta=theta),
         simplify=FALSE,USE.NAMES=TRUE)
}

Zprecompute.lrs.locus <- function(locus,ki.type,fr,theta=0){
  # ladder length
  L <- length(fr[[locus]])
  # all possible geno's
  # make combs (1,1),(2,1),..,(10,1),(2,2),(3,2),..,(10,2),..,(10,10)
  G <- cbind(unlist(sapply(1:L,function(l) l:L)),rep(1:L,L:1))
  colnames(G) <- paste(locus,c(".1",".2"),sep="")
  as.vector(apply(G,1,function(g0) (ki.db(g0,G,hyp.1=ki.type,freqs=fr,theta=theta))))  
}