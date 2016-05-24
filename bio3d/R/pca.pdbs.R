"pca.pdbs" <- function(pdbs, core.find=FALSE, fit=FALSE, ...) {
  ## Log the call
  cl <- match.call()

  if(core.find & fit) {
    warning("incompatible arguments- neglecting 'fit=TRUE'")
    fit=FALSE
  }
  
  if(core.find) {
    core <- core.find(pdbs)
    pdbs$xyz = pdbfit(pdbs, core$c0.5A.xyz)
  } else if(fit) {
     pdbs$xyz = pdbfit(pdbs)
  }
  
  gaps.pos <- gap.inspect(pdbs$xyz)
  pc <- pca.xyz(pdbs$xyz[,gaps.pos$f.inds], ...)

  pc$call=cl
  return(pc)
}
