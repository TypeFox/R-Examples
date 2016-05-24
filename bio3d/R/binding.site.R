"binding.site" <- function(a, b = NULL, a.inds = NULL, b.inds = NULL,
                           cutoff=5, hydrogens=TRUE, byres=TRUE, verbose=FALSE) {

  cl <- match.call()
  sep <- "_"

  trim <- function(s, leading=TRUE, trailing=TRUE) {
    if(leading)
      s <- sub("^ +", "", s)
    if(trailing)
      s <- sub(" +$", "", s)
    s[(s=="")]<-""
    s
  }

  if (!is.pdb(a))
    stop("must supply an input 'pdb' object 'a', i.e. from 'read.pdb'")

  ## workaround for NA chains
  if(any(is.na(a$atom$chain)))
    a$atom$chain[is.na(a$atom$chain)] <- " "

  ## backup of the original pdb provided
  a.orig <- a

  ## two PDBs provided
  if(!is.null(b)) {
    if(!is.pdb(b))
      stop("'b' should be a 'pdb' object as obtained from 'read.pdb'")

    if ( hydrogens ) {
      if(is.null(a.inds))
        a.inds <- atom.select(a, "all", verbose=verbose)
      if(is.null(b.inds))
        b.inds <- atom.select(b, "all", verbose=verbose)
    }
    else {
      if(is.null(a.inds))
        a.inds <- atom.select(a, string='noh', verbose=verbose)
      if(is.null(b.inds))
        b.inds <- atom.select(b, string='noh', verbose=verbose)
    }
  }

  ## one PDB object is provided
  else {
    if(is.null(a.inds) & is.null(b.inds)) {
      a.inds <- atom.select(a, "protein", verbose=verbose)
      b.inds <- atom.select(a, "ligand", verbose=verbose)

      if(!length(a.inds$atom)>0)
        stop("insufficent 'protein' atoms in structure")
      if(!length(b.inds$atom)>0)
        stop("insufficent 'ligand' atoms in structure")
    }

    b <- trim.pdb(a, b.inds)
    a <- trim.pdb(a, a.inds)

    if ( hydrogens ) {
      a.inds <- atom.select(a, "all", verbose=verbose)
      b.inds <- atom.select(b, "all", verbose=verbose)
    }
    else {
      a.inds <- atom.select(a, string='noh', verbose=verbose)
      b.inds <- atom.select(b, string='noh', verbose=verbose)
    }
  }

  if(!(length(a.inds$atom)>0 | length(b.inds$atom)>0))
    stop("insufficent atoms in selection(s)")

  ## omit hydrogens if any
  a <- trim.pdb(a, a.inds)
  b <- trim.pdb(b, b.inds)

  ## Calcuate pair-wise distances
  dmat <- dist.xyz(matrix(a$xyz, ncol=3, byrow=TRUE), matrix(b$xyz, ncol=3, byrow=TRUE))

  ## atoms of a in contact with b
  cmap <- apply(dmat, 1, function(x) any(x <= cutoff))
  atom.inds <- which(cmap)

  ## return NULL if no atoms are closer than cutoff
  if(length(atom.inds)<1)  {
    cat("  no atoms found within", cutoff, "A\n")
    return(NULL)
  }

  ## get rid of any trailing and leading spaces
  a$atom$resid <- trim(a$atom$resid)
  a$atom$resno <- trim(a$atom$resno)
  a$atom$elety <- trim(a$atom$elety)
  a.orig$atom$resno <- trim(a.orig$atom$resno)
  a.orig$atom$elety <- trim(a.orig$atom$elety)

  ## return all atoms in a contacting residue, otherwise, just the atoms
  if(byres) {
    resno.map  <- apply(a$atom[atom.inds, c("resno", "chain")], 1, paste, collapse=sep)
    all.resno  <- apply(a.orig$atom[, c("resno", "chain")],     1, paste, collapse=sep)
    atom.inds2 <- which(all.resno %in% resno.map)
  }
  else {
    resno.map  <- apply(a$atom[atom.inds, c("elety", "resno", "chain")], 1, paste, collapse=sep)
    all.resno  <- apply(a.orig$atom[, c("elety", "resno", "chain")],     1, paste, collapse=sep)
    atom.inds2 <- which(all.resno %in% resno.map)
  }

  xyz.inds <- atom2xyz(atom.inds2)

  ## check for chain IDs
  tmp <- unique(paste(a$atom[atom.inds, "resid"],
                      a$atom[atom.inds, "resno"],
                      a$atom[atom.inds, "chain"],
                      sep=sep))

  resno <- as.numeric(unlist(lapply(strsplit(tmp, sep), function(x) x[2])))
  chain <- unlist(lapply(strsplit(tmp, sep), function(x) x[3]))
  chain[chain==" "] <- NA

  if(all(is.na(chain))) {
    resnames <- unique(paste(a$atom[atom.inds, "resid"], " ",
                             a$atom[atom.inds, "resno"],
                             sep=""))
  }
  else {
    resnames <- unique(paste(a$atom[atom.inds, "resid"], " ",
                             a$atom[atom.inds, "resno"],
                             " (", a$atom[atom.inds, "chain"], ")",
                             sep=""))
  }

  sele <- list(atom=atom.inds2, xyz=xyz.inds)
  class(sele) <- "select"

  out <- list(inds=sele, resnames=resnames, resno=resno, chain=chain, call=cl)
  return(out)
}

