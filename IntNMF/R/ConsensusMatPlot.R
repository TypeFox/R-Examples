ConsensusMatPlot <-
function(fit,rowLab=TRUE,colLab=TRUE){
     if(rowLab==TRUE) row.name <- NULL  else row.name <- NA
     if(colLab==TRUE) col.name <- NULL  else col.name <- NA
     aheatmap(fit$consensus[names(sort(fit$clusters)),names(sort(fit$clusters))],color="Greys",Rowv=NA, Colv=NA,labRow=row.name,labCol=col.name)
}
