stratifyvegdata<-function(x,sizes, plotColumn="plot", speciesColumn = "species", abundanceColumn="abundance", sizeColumn = "size", counts=FALSE, mergeSpecies=FALSE) {
  treeData = as.data.frame(x)
  plotColumnId = which(names(treeData)==plotColumn)
  abundanceColumnId = which(names(treeData)==abundanceColumn)
  sizeColumnId = which(names(treeData)==sizeColumn)
  speciesColumnId = which(names(treeData)==speciesColumn)
  if(mergeSpecies) treeData[,speciesColumnId] = "allspecies"
  spnames =unique(treeData[,speciesColumnId])
 
  stratify<-function(treeDataPlot, sizes, spnames=NULL, speciesColumnId, abundanceColumnId, sizeColumnId, counts=FALSE) {
    if(is.null(spnames)) spnames = unique(treeData[,speciesColumnId])
    nsp = length(spnames)
    nstrata = length(sizes)
    m = data.frame(matrix(0,nrow=nsp, ncol=nstrata))
    row.names(m) = spnames
    if(!is.null(names(sizes))) names(m) = names(sizes)
    else names(m) = paste("S",1:nstrata, sep="")
    for(i in 1:nrow(treeDataPlot)) {
      isp = which(spnames==treeDataPlot[i,speciesColumnId])
      sel = (sizes <= treeDataPlot[i, sizeColumnId])
      if(sum(sel,na.rm=T)>0) {
        istratum = max(which(sel))
        if(!counts) m[isp,istratum] = m[isp,istratum]+treeDataPlot[i,abundanceColumnId]
        else m[isp,istratum] = m[isp,istratum]+1
      }
    }
    return(m)
  }
  
  X = lapply(split(treeData,treeData[,plotColumnId]),
         FUN=stratify,sizes=sizes, 
         spnames=spnames, speciesColumnId = speciesColumnId, abundanceColumnId =abundanceColumnId, sizeColumnId=sizeColumnId, counts=counts)
  class(X)<-c("list","stratifiedvegdata")
  return(X)
}