`pdbseq` <-
function(pdb, inds=NULL, aa1=TRUE) {
  ## b.inds <- atom.select(pdb, "//B////CA/")
  ## seq.pdb(pdb, b.inds)

  if(is.null(inds))
    inds <- atom.select(pdb, "calpha", verbose=FALSE)
#    inds <- which(pdb$calpha)
#    inds <- atom.select(pdb, "//////CA/", verbose=FALSE)$atom

  if(is.list(inds))
    inds <- inds$atom

  if(aa1) {
    aa <- aa321(pdb$atom[inds,"resid"])
  } else {
    aa <- pdb$atom[inds,"resid"]
  }
  if(length(aa) > 0) {
    names(aa) <- pdb$atom[inds,"resno"]
  }
  return(aa)
}

