alleleRto1 <-
function ( geno , marker.label=NULL , miss.val=NA)  {
  geno <- as.matrix(geno)
  if ( !all(is.na(miss.val)) ) {  geno [  geno %in% miss.val ]  <- NA  }
  geno [ !(geno %in% c(1,3,2)) ] <- NA
  geno <- replace(geno,(geno %in% 1),0)
  geno <- replace(geno,(geno %in% 3),1)
  if ( is.null(marker.label) ) marker.label <- paste("S",1:ncol(geno),sep="")
  colnames(geno) <- marker.label
  rownames (geno) <- rownames(geno)
  return ( geno )
}
