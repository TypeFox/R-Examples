#' Counts IBS alleles between case profile and all database profiles
#' 
#' @param x case profile
#' @param db database of profiles
#' @details A profile consists of two alleles at each locus. When two profiles are compared at a locus, there can be 0, 1 or 2 alleles IBS (identical by state). The function compares a single profile with a database of \eqn{N} profiles and counts the number of IBS alleles.
#' @return Data frame with three columns \enumerate{
#' \item \code{ibs}: the number of alleles IBS of \code{x} and a \code{db} profile, ranges from \eqn{0} to \eqn{2*nloci}
#' \item \code{full.matches}: the number of loci for which 2 alleles are IBS, ranges from \eqn{0} to \eqn{nloci}
#' \item \code{partial.matches}: the number of loci for which 1 allele is IBS, ranges from \eqn{0} to \eqn{nloci}
#' }
#' @seealso \code{\link{ibs.pairs}} for comparing many profile pairs
#' @examples
#' ## monte-carlo simulation of number of IBS alleles for FS, PO and UN w.r.t a single case profile
#' 
#' data(freqsNLsgmplus)
#' 
#' #start with a single case profile
#' x <- sample.profiles(N=1,freqs=freqsNLsgmplus)
#' 
#' # sample siblings
#' x.fs <- sample.relatives(x,N=10^3,type="FS") 
#' # sample parent/offspring
#' x.po <- sample.relatives(x,N=10^3,type="PO") 
#' # sample unrelated profiles
#' x.unr <- sample.profiles(N=10^3,freqs=freqsNLsgmplus) 
#' 
#' # make histograms of the number of ibs alleles
#' hist(ibs.db(x,x.fs)$ibs,xlim=c(0,20),main="FS",xlab="IBS")
#' hist(ibs.db(x,x.po)$ibs,xlim=c(0,20),main="PO",xlab="IBS")
#' hist(ibs.db(x,x.unr)$ibs,xlim=c(0,20),main="UN",xlab="IBS")
#'
#' @export
ibs.db <- function(x,db){
  x <- Zassure.matrix(x)
  
  #some error checking
  if (nrow(x)>1) warning("nrow(x)>1, only first profile is used!")
  x.loci <- Znames.to.loci(Zprofile.names(x))
  db.loci <- Znames.to.loci(colnames(db))
  if (!all(x.loci %in% db.loci)) warning("not all loci of target profile are contained in db")
  
  #loop through loci and count the ibs alleles
  ibs.2 <- ibs.1 <- ibs.2.loc <- ibs.1.loc <- rep(0L,nrow(db))
  for (locus.i in seq(x.loci)){
    locus.name <- x.loci[locus.i]
    
    if (locus.name %in% db.loci){
      ind <- locus.i*2+c(-1,0)
      
      a <- x[1,ind[1]];      b <- x[1,ind[2]]
      c <- db[,paste(locus.name,".1",sep="")] #db
      d <- db[,paste(locus.name,".2",sep="")]
      
      ac <- (a==c);    ad <- (a==d)
      bc <- (b==c);    bd <- (b==d)
      
      ac[is.na(ac)] <- FALSE
      ad[is.na(ad)] <- FALSE
      bc[is.na(bc)] <- FALSE
      bd[is.na(bd)] <- FALSE
      
      ibs.2.loc <- as.integer((ad&bc)|(ac&bd))
      ibs.1.loc <- as.integer(xor(ac,bd)|xor(ad,bc))    
      
      ibs.2 <- ibs.2 + ibs.2.loc 
      ibs.1 <- ibs.1 + ibs.1.loc
    }
  }
  data.frame(ibs=(2L*ibs.2+ibs.1),full.matches=ibs.2,partial.matches=ibs.1)
}
NULL
#' Counts IBS alleles between profile pairs
#' 
#' @param x1 N profiles
#' @param x2 N profiles
#' @details A profile consists of two alleles at each locus. When two profiles are compared at a locus, there can be 0, 1 or 2 alleles IBS (identical by state). The function compares two databases of profiles (\code{x1} and \code{x2}) and counts the number of IBS alleles between profile pairs. The function expects \code{x1} and \code{x2} to be of equal size.
#' @return Data frame with three columns, containing for each profile pair:\enumerate{
#' \item \code{ibs}: the number of alleles IBS, ranges from \eqn{0} to \eqn{2*nloci}
#' \item \code{full.matches}: the number of loci for which 2 alleles are IBS, ranges from \eqn{0} to \eqn{nloci}
#' \item \code{partial.matches}: the number of loci for which 1 allele is IBS, ranges from \eqn{0} to \eqn{nloci}
#' }
#' @seealso \code{\link{ibs.db}} for comparing one profile against a \code{db} 
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
#' legend("topright",legend=c("PO","UN"),fill=c("red","blue"))
#' @export
ibs.pairs <- function(x1,x2){
  #some error checking
  if (nrow(x1)!=nrow(x2)) stop("x1 and x2 do not contain the same number of profiles")
  x1.loci <- Znames.to.loci(colnames(x1))
  x2.loci <- Znames.to.loci(colnames(x2))
  if (!all(x1.loci %in% x2.loci)) warning("not all loci of x2 profile are contained in x2")
  
  ibs.2 <- ibs.1 <- ibs.2.loc <- ibs.1.loc <- rep(0L,nrow(x1))
  for (locus.i in 1:(ncol(x1)/2)){
    locus.name <- x1.loci[locus.i]
    if (locus.name %in% x2.loci){    
      ind <- locus.i*2+c(-1,0)
      a <- x1[,ind[1]];    b <- x1[,ind[2]]
      c <- x2[,paste(locus.name,".1",sep="")]
      d <- x2[,paste(locus.name,".2",sep="")]
      
      ac <- (a==c);    ad <- (a==d)
      bc <- (b==c);    bd <- (b==d)
      
      ac[is.na(ac)] <- FALSE
      ad[is.na(ad)] <- FALSE
      bc[is.na(bc)] <- FALSE
      bd[is.na(bd)] <- FALSE
      
      ibs.2.loc <- as.integer((ad&bc)|(ac&bd))
      ibs.1.loc <- as.integer(xor(ac,bd)|xor(ad,bc))    
  
      ibs.2 <- ibs.2 + ibs.2.loc 
      ibs.1 <- ibs.1 + ibs.1.loc
    }
  }
  data.frame(ibs=(2L*ibs.2+ibs.1),full.matches=ibs.2,partial.matches=ibs.1)  
}
NULL