## Use for trimming a pdbs object, either by removing structures,
## or by removing columns 

trim.pdbs <- function(pdbs, row.inds=NULL, col.inds=NULL, ...) {
  if(!inherits(pdbs, "pdbs"))
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")
  
  ## Log the call
  cl <- match.call()
  
  if(is.null(row.inds))
    row.inds <- seq(1, nrow(pdbs$resno), by=1)
  if(is.null(col.inds)) {
    gaps <- gap.inspect(pdbs$resno[row.inds,,drop=FALSE])
    col.inds <- which(gaps$col < dim(pdbs$resno[row.inds,,drop=FALSE])[1L])
  }
  
  if(any(col.inds<0))
    col.inds.xyz <- atom2xyz(abs(col.inds)) * sign(rep(col.inds, each=3))
  else
    col.inds.xyz <- atom2xyz(col.inds)
  
  new <- NULL
  new$id    =pdbs$id[row.inds]
  new$xyz   =pdbs$xyz[row.inds, col.inds.xyz, drop=FALSE]
  new$resno =pdbs$resno[row.inds, col.inds, drop=FALSE]
  new$b     =pdbs$b[row.inds, col.inds, drop=FALSE]
  new$chain =pdbs$chain[row.inds, col.inds, drop=FALSE]
  new$ali   =pdbs$ali[row.inds, col.inds, drop=FALSE]
  new$resid =pdbs$resid[row.inds, col.inds, drop=FALSE]
  new$sse   =pdbs$sse[row.inds, col.inds, drop=FALSE]
  new$call  =cl
 
  class(new) <- c("pdbs", "fasta")
  return(new)
}
