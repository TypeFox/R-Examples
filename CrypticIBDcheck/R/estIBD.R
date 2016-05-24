IBDest.study <- function(snpobjtr, cdlibs) {
   #remove monomorphic SNPs, which yield missing values in all of the 
   # non-zero cdl IBS probabilities, such as P(I=0 | Z=0) in the first
   # column of cldibs.
   excl<-is.na(cdlibs[,1]) 
   snpobjtr<-snpobjtr[,!excl]
   cdlibs<-cdlibs[!excl,]
   out <- .Call("IBDest_study", t(snpobjtr@.Data), 
	nrow(snpobjtr),
	ncol(snpobjtr), 
	t(cdlibs),
	new.env(), PACKAGE="CrypticIBDcheck")
   pz0 <- matrix(out[[1]], nrow(snpobjtr), nrow(snpobjtr), byrow=TRUE)
   pz1 <- matrix(out[[2]], nrow(snpobjtr), nrow(snpobjtr), byrow=TRUE)
   pz2 <- matrix(out[[3]], nrow(snpobjtr), nrow(snpobjtr), byrow=TRUE)
   # following commands vec-out the upper-triangular matrix column 
   # by column.
   pz0 <- pz0[upper.tri(pz0)]
   pz1 <- pz1[upper.tri(pz1)]
   pz2 <- pz2[upper.tri(pz2)]
   data.frame(pz0=pz0, pz1=pz1, pz2=pz2)
}

IBDest.sim <- function(snpmat, cdlibs) {
   if(is.null(snpmat)) {
     return(NULL)
   }
   excl<-is.na(cdlibs[,1]) 
   snpmat<-snpmat[,!excl]
   cdlibs<-cdlibs[!excl,]
   out <- .Call("IBDest_sim", t(snpmat@.Data), 
	nrow(snpmat)/2,
	ncol(snpmat), 
	t(cdlibs),
	new.env(), PACKAGE="CrypticIBDcheck")
  names(out) <- c("pz0", "pz1", "pz2")
  data.frame(out)
}
