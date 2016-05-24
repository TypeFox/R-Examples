"atom.select.prmtop" <- function(prmtop, ...) {
  if(!inherits(prmtop, "prmtop"))
    stop("provide a PRMTOP object as obtained from read.prmtop()")

  cl <- match.call()
  tmp.pdb <- as.pdb.prmtop(prmtop, crd=as.numeric(rep(NA, prmtop$POINTERS[1]*3)))
  sele <- atom.select.pdb(tmp.pdb, ...)
  sele$call <- cl
  return(sele)
}
