CheckPrintGeneListFormat <- function(lpc.obj, numGenes=100, gene.names=NULL, lpcfdr.out=NULL){
  if(numGenes<1 || floor(numGenes)!=numGenes) stop("numGenes must be an integer.....")
  if(length(lpc.obj$lpcscores)!=length(lpc.obj$tscores)) stop("LPC and T scores must be the same length...")
  if(!is.null(lpcfdr.out)){
    if(length(lpcfdr.out$fdrt)!=length(lpcfdr.out$fdrlpc) || length(lpcfdr.out$fdrt)!=length(lpc.obj$lpcscores)) stop("There must be one LPC FDR and one T FDR for each gene")
  }
  if(!is.null(gene.names) && length(gene.names)!=length(lpc.obj$lpcscores)) stop("There must be the same number of gene names as there are gene scores.")
}

