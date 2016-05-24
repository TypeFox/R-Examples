layout.cna <- function(x, pdb, renumber=TRUE, k=2, full=FALSE){

  ## Return the coordinate centers of network communities 
  ##  as defined in "x$communities$membership" 'membership vector' 
  ##  using Calpha's in 'pdb'.
  ##
  ##  If k=3 the xyz geometric centers are returned. if k<3 then
  ##  multidimensional scaling is used for k space ordination.
  ##
  ##   co2  <- layout.cna(net, pdb)
  ##   co3  <- layout.cna(net, pdb, k=3)
  ##   all2 <- layout.cna(net, pdb, k=2, full=TRUE)
  ##   plot.cna(net, layout=co2)

  if( !inherits(pdb, "pdb") ) {
    stop("Input 'pdb' is not of class 'pdb' as obtained from 'read.pdb()'")
  }
  if(!k %in% c(1,2,3)) {
    stop("Input 'k' should have a value of 3, 2 or 1")
  }
  if(inherits(x, "cna") || is.list(x)) {
    if(!full) {
      ## We want community coords
      membership <- x$communities$membership
    } else {  
      ## We want full all-atom/Calpha coords (check network for number)
      membership <- c(1:length(x$communities$membership))
    }
  } else {
    if(full) {
      ## Assuming we want Calpha coords - this is a BIG assumption! 
      membership <- c(1:sum(pdb$calpha))
    }
    stop("Input object 'x' should be of class cna")
  }

  ## Renumber 'pdb' to match membership resno indices
  if(renumber) {
    pdb <- convert.pdb(pdb, renumber=TRUE, rm.h=FALSE, verbose=FALSE)
  }

  ##-- Check if the number of residues in 'pdb' equals
  ##   the length of 'membership' vector
  notprotein.inds <- atom.select(pdb, "notprotein", verbose=FALSE)

  if(length(notprotein.inds$atom)>0){
    num.res <- length(pdb$atom[pdb$calpha,"resno"]) + length(unique(pdb$atom[notprotein.inds$atom,7]))
  }
  if(length(notprotein.inds$atom)==0){
    num.res <- length(pdb$atom[pdb$calpha,"resno"])
  }
  
  ##-- Calculate the geometric center of each community
  n <- unique(membership[!is.na(membership)])
  
  cent <- matrix(NA, nrow=length(n), ncol=3)
  a <- 1
  for(i in n){
    inds <- atom.select(pdb, resno=which(membership==i),
                        elety="CA", verbose=FALSE)

    cent[a,] <- apply( matrix(pdb$xyz[inds$xyz], nrow=3), 1, mean)
    a <- a + 1 
  }

  if(k != 3) {
    ##-- Multidimensional scaling for 2D or 1D projection
    ##   note. dist(centers) and dist.xyz(centers) give same answer
    if(nrow(cent) - 1 < k) {
       # e.g. only two communities 
       cent <- cmdscale(dist(cent), k = nrow(cent) - 1)
       cent <- cbind(cent, matrix(0, nrow=nrow(cent), ncol=k-nrow(cent)+1))
    } else { 
       cent <- cmdscale(dist(cent),k=k)
    }
  }
  return(cent)
}
