dssp.xyz <- function(xyz, pdb, ...) {
  if(!is.pdb(pdb))
    stop("provide a pdb object as obtained from function 'read.pdb'")

  if(!is.xyz(xyz) && !is.matrix(xyz))
    stop("provide an xyz object containing the trajectory coordinates")

  sse.mat <- NULL
  dims <- dim(xyz)
  for (i in 1:dims[1L]) {
    pdb.tmp     <- pdb
    pdb.tmp$xyz <- as.xyz(xyz[i,])
    sse     <- dssp.pdb(pdb.tmp, ...)$sse
    sse.mat <- rbind(sse.mat, sse)
  }

  ##sse.mat[ sse.mat==" " ] <- "-"
  return(sse.mat)
}
