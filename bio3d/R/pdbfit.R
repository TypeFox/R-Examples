pdbfit <- function(...)
  UseMethod("pdbfit")

pdbfit.pdb <- function(pdb, inds=NULL, ...) {
  if(!is.pdb(pdb))
    stop("Input 'pdb' should be of class 'pdb', e.g. from read.pdb()")

  if(nrow(pdb$xyz)<2)
    stop("nothing to fit. < 2 frames in pdb object")

  if(is.null(inds)) {
    inds <- atom.select(pdb, "calpha")
    cat("  no indices provided. using all", length(inds$atom), "calpha atoms\n")
  }

  if(length(inds$xyz)<3)
    stop("insufficent atoms to superimpose")

  return(fit.xyz( fixed=pdb$xyz[1,], mobile=pdb$xyz,
                 fixed.inds=inds$xyz, mobile.inds=inds$xyz))
  
}

pdbfit.pdbs <- function(pdbs, inds=NULL, outpath=NULL, ...) {
  ##
  ## Quick Fit Fitter for PDBs
  ##  was called 'fit.pdbs()' in model.R
  ##
  if(!inherits(pdbs, "pdbs")) {
    stop("Input 'pdbs' should be of class 'pdbs', e.g. from pdbaln() or read.fasta.pdb()")
  }
  full <- ifelse(is.null(outpath), FALSE, TRUE)
  if(is.null(inds)) {    inds <-gap.inspect(pdbs$xyz)$f.inds  }
  if(is.list(inds)){ inds=inds$xyz }
  invisible( fit.xyz( fixed=pdbs$xyz[1,], mobile=pdbs, fixed.inds=inds,
                  mobile.inds=inds, outpath=outpath,
                  full.pdbs=full, ... ))
}
