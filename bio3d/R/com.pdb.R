"com.pdb" <-
  function(pdb, inds=NULL, use.mass=TRUE, ... ) {
    
    if (missing(pdb))
      stop("Please supply an input 'pdb' object, i.e. from 'read.pdb()'")
    if(!is.pdb(pdb))
      stop("Input 'pdb' must be of type 'pdb'")
    
    if(is.null(inds)) {
      xyz <- pdb$xyz
      at <- pdb$atom[, "elety"]
    }
    else {
      if(!is.select(inds))
        stop("provide a select object as obtained from 'atom.select'")
      
      if(length(inds$xyz)<3)
        stop("insufficient atoms in selection")
      xyz <- pdb$xyz[,inds$xyz]
      at <- pdb$atom[inds$atom, "elety"]
    }
    
    if(use.mass) {
      m <- atom2mass(at, ...)
    }
    else {
      m <- NULL
    }

    com <- com.xyz(xyz, m)
    return(com)
  }
