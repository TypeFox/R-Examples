as.pdb.prmtop <- function(prmtop, crd=NULL, inds=NULL, inds.crd=inds, ncore=NULL, ...) {
  ncore <- setup.ncore(ncore, bigmem=FALSE)

  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply

  if(!inherits(prmtop, "prmtop"))
    stop("provide a prmtop object as obtained from function read.prmtop")

  natoms.prmtop <- prmtop$POINTERS[1]
  if(is.null(crd)) {
    warning("producing PDB object with no XYZ coordinates")
    crd <- as.numeric(rep(NA, natoms.prmtop*3))
  }

  if(!inherits(crd, "amber")) {
    new <- list()
    new$xyz <- as.xyz(crd)
    new$natoms <- ncol(new$xyz)/3
    crd <- new
  }

  natoms.crd <- crd$natoms
  if( any(c(!is.null(inds), !is.null(inds.crd))) ) {
    if(is.null(inds)) {
      inds$atom = seq(1, natoms.prmtop)
      inds$xyz = atom2xyz(inds$atom)
      class(inds) = "select"
    }

    if(is.null(inds.crd)) {
      inds.crd$atom = seq(1, natoms.crd)
      inds.crd$xyz = atom2xyz(inds.crd$atom)
      class(inds.crd)="select"
    }

    natoms.prmtop = length(inds$atom)
    natoms.crd    = length(inds.crd$atom)
  }

  if(natoms.prmtop != natoms.crd)
    stop(paste("atom number mismatch:", natoms.prmtop, "(prmtop) vs", natoms.crd, "(crds)"))

  resmap <- function(i, type='resid') {
    if(i==length(prmtop$RESIDUE_POINTER))
      j <- prmtop$POINTERS[1] - prmtop$RESIDUE_POINTER[i] + 1
    else
      j <- prmtop$RESIDUE_POINTER[i+1] - prmtop$RESIDUE_POINTER[i]
    if(type=='resno')
      return(rep(i,j))
    if(type=='resid')
      return(rep(prmtop$RESIDUE_LABEL[i], j))
  }

  resno <- unlist(mylapply(1:length(prmtop$RESIDUE_POINTER), resmap, 'resno'))
  resid <- unlist(mylapply(1:length(prmtop$RESIDUE_POINTER), resmap, 'resid'))

  if(any(c(!is.null(inds), !is.null(inds.crd)))) {
    pdb <- as.pdb.default(xyz=crd$xyz[,inds.crd$xyz], elety=prmtop$ATOM_NAME[inds$atom],
                          resno=resno[inds$atom], chain=as.character(NA), resid=resid[inds$atom])
  }
  else {
    pdb <- as.pdb.default(xyz=crd$xyz, elety=prmtop$ATOM_NAME,
                          resno=resno, chain=as.character(NA), resid=resid)
  }

  pdb$call = match.call()
  return(pdb)
}
