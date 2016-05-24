# Automatically generated from all.nw using noweb
 
## renamed from pedBits, part of pedigree.shrink functions

bitSize <- function(ped) {
  ## calculate bit size of a pedigree
  
  if(class(ped) != "pedigree")
    stop("Must be a pegigree object.\n")
  
  father = ped$findex
  mother = ped$mindex
  id = ped$id
  
  founder <- father==0 & mother==0
  pedSize <- length(father)
  nFounder <- sum(founder)
  nNonFounder <- pedSize - nFounder
  bitSize <- 2*nNonFounder - nFounder
  return(list(bitSize=bitSize,
              nFounder = nFounder,
              nNonFounder = nNonFounder))
}

