cmap.pdb <- function(pdb, inds=NULL, verbose=FALSE, ...) {
  if(!is.pdb(pdb))
    stop("provide a pdb object as obtained from function 'pdb'")

  if(is.null(inds)) {
    inds <- atom.select(pdb, "notwater", verbose=verbose)
  }

  pdb <- trim.pdb(pdb, inds)
  xyz <- pdb$xyz
  grpby <- paste(pdb$atom$chain, pdb$atom$insert, pdb$atom$resno, sep="-")
  return(cmap.xyz(xyz, grpby, ...))
  

}
