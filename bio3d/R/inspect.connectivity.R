## Useful for checking the connectivity in a pdb(s) object

"inspect.connectivity" <- function(pdbs, cut=4.) {
  xyz <- NULL; ids <- NULL;
  if(inherits(pdbs, "pdbs")) {
    xyz <- pdbs$xyz
    n <- length(pdbs$id)
    ids <- pdbs$id
  }
  else if(is.pdb(pdbs)) {
    ca.inds <- atom.select(pdbs, 'calpha', verbose=FALSE)
    xyz <- as.xyz(pdbs$xyz)[1, ca.inds$xyz, drop=FALSE]
    n <- 1
  }
  else if(inherits(pdbs, "xyz")) {
    xyz <- pdbs
    n <- nrow(xyz)
  }
  else {
    stop("Please provide coordinates as a \n 'pdbs', 'pdb', or xyz matrix format")
  }

  if(length(xyz)<6) {
    warning("Insufficient C-alpha atoms in structure to determine connectivity")
    return(FALSE)
  }

  is.connected <- function(xyz) {
    xyz <- matrix(xyz[!is.na(xyz)], ncol=3, byrow=T)
    for(i in 1:(nrow(xyz)-1)) {
      d <- sqrt((xyz[i,1]-xyz[i+1,1])**2 +
                (xyz[i,2]-xyz[i+1,2])**2 +
                (xyz[i,3]-xyz[i+1,3])**2 )

      if(d>cut)
        return(FALSE)
    }
    return(TRUE)
  }

  cons <- rep(NA, length=n)
  for(i in 1:n) {
    cons[i] <- is.connected(xyz[i,])
  }

  names(cons) <- ids
  return(cons)
}
