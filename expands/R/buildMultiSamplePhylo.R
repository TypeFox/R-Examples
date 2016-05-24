buildMultiSamplePhylo<-function (samGr,out,treeAlgorithm="bionjs", ambigSg=F,plotF=1,spRes=1){
  library(expands)
  if(!is.na(spRes) && !spRes){
    print("Warning: Calculating cross-sample phylogeny at metapopulation resolution")
  }
  cols = c("Count", "chr", "startpos", "endpos",  "CN_Estimate")
  dummySNVcols=c("Count","endpos");
  allCBS = c()
  allDM = c()
  dmPris=list();
  n_Samples = length(samGr$labels)
  for (i in 1:n_Samples) {
    cbsPri = read.table(samGr$cbs[[i]], sep = "\t", 
                        header = T, stringsAsFactors = FALSE)
    cbsPri[, c("startpos", "endpos")] = round(cbsPri[, c("startpos", 
                                                         "endpos")]/10000) * 10000
    dmPri = read.table(samGr$sps[[i]], sep = "\t", 
                       header = T, stringsAsFactors = FALSE)
    tmpI=grep("SP_0", colnames(dmPri));
    if (!isempty(tmpI)){
        dmPri=dmPri[,-1*tmpI]; ##remove SP composition specific columns
    }
    if(!is.na(spRes) && !spRes){
      dmPri[,"SP"]=max(dmPri[,"SP"],na.rm=T);
      if(any("SP_cnv" %in% colnames(dmPri))){ ##Backward compatibility
        dmPri[,"SP_cnv"]=max(dmPri[,"SP"],na.rm=T);
      }
    }
    for (j in 1:length( dummySNVcols)){
      if (!any(colnames(dmPri)==dummySNVcols[j])){
        tmp=colnames(dmPri);
        dmPri=cbind(dmPri,matrix(NA,nrow(dmPri),1));
        colnames(dmPri)=c(tmp,dummySNVcols[j]);
        if(dummySNVcols[j]=='endpos'){
          dmPri[,dummySNVcols[j]]=dmPri[,'startpos'];
        }
      }
    }
    
    ##Add optional columns if they don't exist
    cbsPri = .addMissingCols(cbsPri)
    dmPri = .addMissingCols(dmPri)
    dmPris[[i]] = dmPri
    allCBS = as.matrix(rbind(allCBS, cbsPri))
    allDM = as.matrix(rbind(allDM, dmPri))
    allCBS[, "Count"] = c(1:nrow(allCBS))
    allDM[, "Count"] = c(1:nrow(allDM))
  }
  dupI = which(duplicated(allCBS[, c("chr", "startpos", "endpos")]))
  if (length(dupI) > 0) {
    allCBS = allCBS[-1 * dupI, ]
  }
  dupI = which(duplicated(allDM[, c("chr", "startpos")]))
  if (length(dupI) > 0) {
    allDM = allDM[-1 * dupI, ]
  }
  aqCBS = allCBS[, cols]
  aqDM = allDM[, cols]
  for (i in 1:n_Samples) {
    dmPri = dmPris[[i]];
    print(paste("Processing sample ", i, " out of ",n_Samples,sep=""));
    aQpriCBS = try(assignQuantityToSP(allCBS[, cols], dmPri,  colName = "PM",
                                      keepAmbigSeg = ambigSg), silent = FALSE)
    aQpriCBS=aQpriCBS$ploidy;
    dmPri[, "PM_B"] = sign(dmPri[, "PM_B"])
    aQpriDM = try(assignQuantityToSP(allDM[, cols], dmPri, 
                                     colName = "PM_B"), silent = FALSE)
    aQpriDM=aQpriDM$ploidy;
    firstI = min(grep("SP", colnames(aQpriCBS)))
    aqCBS = cbind(aqCBS, aQpriCBS[, firstI:ncol(aQpriCBS)])
    aqDM = cbind(aqDM, aQpriDM[, firstI:ncol(aQpriDM)])
    nSPs = length(unique(dmPri[!is.na(dmPri[, "SP"]), "SP"]))
    lab = paste(samGr$labels[[i]], "_SP", sep = "")
    colns = colnames(aQpriDM)
    colnames(aqCBS) = c(colnames(aqCBS[, 1:(ncol(aqCBS) - 
                                              nSPs)]), gsub("SP", lab, colns[firstI:ncol(aQpriDM)]))
  }
  aQ = rbind(aqCBS, aqDM)
  tr = NULL
  if (class(aQ) == "try-error" || is.null(ncol(aQ))) {
    print("Error encountered while reconstructing phylogeny")
  }
  else {
    trout = buildPhylo(aQ, out, treeAlgorithm = treeAlgorithm)
    tr=trout$tree;
    if (plotF>0){
      jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      colmap = jet(n_Samples)
      colors <- rep(colmap[1], each = length(tr$tip.label))
      for (i in 1:n_Samples) {
        ii = grep(samGr$labels[[i]], tr$tip.label)
        colors[ii] = colmap[i]
      }
      plot(tr, tip.col = colors, cex = 1, type = "u")
    }
    return(tr)
  }
  return(NULL)
}

.addMissingCols<-function(tmp){
  if(!"Count" %in% colnames(tmp)){
    cols=colnames(tmp);
    tmp=cbind(tmp,matrix(NA,nrow(tmp),1))
    colnames(tmp)=c(cols,"Count");
  }
  if(!"SP_cnv" %in% colnames(tmp) && "SP" %in% colnames(tmp)){
    cols=colnames(tmp);
    tmp=cbind(tmp,tmp[,"SP"])
    colnames(tmp)=c(cols,"SP_cnv");
  }
  return(tmp);
}
