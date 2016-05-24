PrintGeneList <- function(lpc.obj, numGenes=100, gene.names=NULL, lpcfdr.out=NULL){
  CheckPrintGeneListFormat(lpc.obj,numGenes,gene.names,lpcfdr.out)
  if(is.null(gene.names)) gene.names <- paste("Gene", sep="", as.character(1:length(lpc.obj$lpcscores))) 
  ord <- order(abs(lpc.obj$lpcscores), decreasing=T)
  if(is.null(lpcfdr.out)){ 
    printframe <- data.frame(cbind(gene.names[ord[1:numGenes]], round(lpc.obj$lpcscores[ord[1:numGenes]],4), round(lpc.obj$tscores[ord[1:numGenes]],4))) 
    names(printframe) <- c("Gene Name", "LPC Score", "T Score")
  } else {  
    printframe <- data.frame(cbind(gene.names[ord[1:numGenes]], round(lpc.obj$lpcscores[ord[1:numGenes]],4), round(lpcfdr.out$fdrlpc[ord[1:numGenes]],4), round(lpc.obj$tscores[ord[1:numGenes]],4), round(lpcfdr.out$fdrt[ord[1:numGenes]],4))) 
    names(printframe) <- c("Gene Name", "LPC Score", "LPC FDR", "T Score", "T FDR") 
  }
  print(printframe) 
}
