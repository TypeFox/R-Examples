"dm" <- function(...)
  UseMethod("dm")

"dm.pdb" <- function(pdb, inds=NULL, grp=TRUE, verbose=TRUE, ...) {
  if(!is.pdb(pdb)) {
    stop("input 'pdb' should be either:
            1. an object returned from from 'read.pdb' or
            2. a numeric 'xyz' vector of coordinates")
  }

  if(!is.null(inds)) {
    pdb <- trim.pdb(pdb, inds)
  }

  if(grp)
    grpby <- paste(pdb$atom$resno, pdb$atom$chain, sep="-")
  else
    grpby <- NULL
  
  d <- dm.xyz(pdb$xyz, grpby=grpby, ...)
  
  class(d) <- "dmat"
  return(d)
}

