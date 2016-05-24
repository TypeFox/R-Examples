runExPANdS<-function(SNV, CBS, maxScore=2.5, max_PM=6, min_CellFreq=0.1, precision=NA, plotF=2,snvF=NULL,maxN=8000,region=NA, peakselection='localsum'){
  if (!exists("SNV") || !exists("CBS")){
    print("Input-parameters SNV and/or CBS missing. Please add the paths to tabdelimited files containing the SNVs and copy numbers.");
    return();
  }
  if (is.na("maxScore")){
    maxScore=100;
  }
  
  if (is.na("min_CellFreq")){
    min_CellFreq=0;
  }
  
  if (is.na("max_PM")){
    print("Parameter max_PM must be set")
    return();
  }
  
  if (is.na("plotF")){
    plotF=0;
  }
  
  nullResult=list("finalSPs"=NULL,"dm"=NULL,"densities"=NULL,"ploidy"=NULL);
  dirF=getwd();
  ##SNVs
  if (is.character(SNV) && file.exists(SNV)){
    print(paste("Running ExPANdS on: ",SNV))
    dm=read.table(SNV,sep="\t",header=TRUE,stringsAsFactors = FALSE);
    if (!all(sapply(dm,is.numeric))){
      print(paste("Warning: not all columns in", SNV,"are numeric. Use only numeric data as input to ensure unexpected conversion does not occur."))
    }
    dm <- dm[ as.character(dm[,"chr"]) %in% as.character(seq(100)), ];
    dm=data.matrix(dm);
    print("Only SNVs with autosomal coordinates included.")
    if(is.null(snvF)){
      snvF=SNV;
    }
  }else if (is.matrix(SNV)){
    dm=SNV;
    if (!is.numeric(dm)) {
      print("SNV matrix has to be numeric. Likely cause: only mutations detected on autosomes accepted for ExPANdS model. Remove SNVs with allosomal and mitochondrial coordinates before you proceed.")
      return();
    }
    if (is.null(snvF)){
      snvF="out.expands";
    }
  }else{
    print("No SNVs provided. Aborting ExPANdS.");
    return();
  }
  
  ##Output
  dirF=fileparts(snvF)$pathstr;
  if (nchar(dirF)==0){
    dirF=".";
  }
  snvF=paste(fileparts(snvF)$name,fileparts(snvF)$ext,sep="");
  
  ##CBS
  if (is.character(CBS) && file.exists(CBS)){
    copyNumber=as.matrix(read.table(CBS,sep="\t",header=TRUE,stringsAsFactors = FALSE))
    if (!all(sapply(copyNumber,is.numeric))){
      print(paste("Warning: not all columns in", CBS,"are numeric. Use only numeric data as input to ensure unexpected conversion does not occur."))
    }
  }else if (is.matrix(CBS)){
    copyNumber=CBS;
  }else{
    print("No copy number information provided. Aborting ExPANdS.")
    return();
  }
  if(!any("CN_Estimate" %in% colnames(dm))){  
    if (any(copyNumber[,"CN_Estimate"]<0) || quantile(copyNumber[,"CN_Estimate"],0.9)<1){
      print("Column <CN_Estimate> in CBS input seems to contain log-ratio entries. Please supply absolute copy number values (e.g. average ~2.0 expected for predominantly diploid genomes).")
      return();
    }
    if (sum(abs(copyNumber[,"CN_Estimate"]-round(copyNumber[,"CN_Estimate"])))<1){
      print("Warning! Copy numbers have values rounded to the closest integer. Using rational positive estimates of copy numbers is recommended.")
    }
    dm=assignQuantityToMutation(dm,copyNumber,"CN_Estimate");
  }else{
    print("Using column <CN_Estimate> from <SNV> as copy number estimate. Parameter <CBS> used only for phylogeny inference.")
  }
  
  ii=which(is.na(dm[,"CN_Estimate"]));
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to unavailable copy number in that region."));
    dm=dm[-ii,];
  }
  ii=which(dm[,"CN_Estimate"]<1);
  homDelRegions=c();
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to homozygous deletions within that region."));
    homDelRegions=dm[ii,];
    dm=dm[-ii,];
  }
  ii=which(dm[,"CN_Estimate"]>max_PM);
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to high-level amplifications (>",max_PM, "copies) within that region. Consider increasing value of parameter max_PM to facilitate inclusion of these SNVs, provided high coverage data (> 150 fold) is available"));
    dm=dm[-ii,];
  }
  ii=which(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]<min_CellFreq);
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to AF*CN below ", min_CellFreq," (SNV can't be explained by an SP present in ",min_CellFreq*100 ,"% or more of the sample)."));
    dm=dm[-ii,];
  }
  
  
  if(nrow(dm)<20){
    print("Not enough mutations provided. Minimum 20 SNVs required to attempt a run.")
    return(nullResult);
  }
  
  ##Regions of interest
  if (is.character(region) && file.exists(region)){
    roi=read.table(region,sep="\t",header=TRUE,stringsAsFactors = FALSE);
    roi <- roi[ as.character(roi[,"chr"]) %in% as.character(seq(100)), ];
    roi=data.matrix(roi);
  }else if (is.matrix(region)){
    roi=region;
  }
  
  idx_R=c(1:size(dm,1));
  if (size(dm,1)>maxN){
    if (is.na(region)){
      print(paste("Input contains more than ",maxN," SNVs and parameter <region> not set. Using default regional boundary (SureSelectExome_hg19)"));
    }
    
    if (!is.numeric(roi)) {
      print("Region matrix has to be numeric.");
      return(nullResult);
    }
    print(paste("Gathering SNVs within regions of interest ..."))
    ##Find SNVs within regions of interest
    idx_R=c();
    so=sort(roi[,'start'],index.return=T);
    roi=roi[so$ix,];
    for(i in 1:nrow(dm)){
      ii=which(roi[,'chr']==dm[i,'chr'])
      ix=which.min(abs(roi[ii,'start']-dm[i,'startpos']))
      ii=ii[c(max(0,(ix - 2)):min(length(ii),(ix + 2)))];
      if(i %% 100 ==0){
        print(paste(i," out of ",nrow(dm), " SNVs tested"))
      }
      for (j in ii){
        if (roi[j,'start']<=dm[i,'startpos'] && roi[j,'end']>=dm[i,'startpos']){
          idx_R=cbind(idx_R,i);
          break;
        }
      }
    }
    print(paste('Found ',length(idx_R), ' SNVs within regions of interest',sep=""))
    idx_R=sample(idx_R,min(maxN,length(idx_R)),replace=FALSE);
    print(paste("Keeping ",length(idx_R), ' of these SNVs (randomly selected).',sep=""))
  }else{
    if (!is.na(region)){
      print(paste("Input contains less than ",maxN," SNVs. Parameter <region> ignored."));
    }
  }
  if (is.na(precision)){
    precision=0.1/log(length(idx_R)/7);
  }
  print(paste('Peak selection strategy: ',peakselection))
  
  #############################
  ###Input parsing ends here###
  #dm=dm[idx_R,]
  cfd=computeCellFrequencyDistributions(dm, max_PM, precision, min_CellFreq=min_CellFreq);
  toUseIdx=which(apply(is.finite(cfd$densities),1,all) );
  toUseIdx=intersect(toUseIdx,idx_R);
  SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision, min_CellFreq=min_CellFreq);
  if(is.null(SPs) || size(SPs,1)==0 || all(is.na(SPs))){
    print("No SPs found.")
    result=list("finalSPs"=NULL,"dm"=cfd$dm,"densities"=cfd$densities);
    return(result);
  }
  if(is.null(dim(SPs)) || nrow(SPs)==1){
    finalSPs=SPs[SPs["score"]<=maxScore];
  }else{
    finalSPs=SPs[SPs[,"score"]<=maxScore,];
  }
  if(!is.null(dim(finalSPs)) && nrow(finalSPs)>1){
    ia=order(finalSPs[,"Mean Weighted"]);
    finalSPs=matrix(finalSPs[ia,],nrow=nrow(finalSPs), ncol=ncol(finalSPs), dimnames=list(1:nrow(finalSPs),colnames(finalSPs)));
  }
  if(size(finalSPs,1)==0){
    print(paste("No SPs found below score:",maxScore))
    result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities);
    return(result);
  }
  ###########################
  ###Assigning SNVs to SPs###
  aM = assignMutations( cfd$dm, finalSPs,peakselection=peakselection);
  dm=aM$dm; SPs=aM$finalSPs;
  
  ###########################
  ##Choose between doublets##
  continueprune=TRUE;
  while (!is.null(dim(SPs)) && nrow(SPs)>1 && continueprune){
    x=sort(SPs[,'Mean Weighted'],index.return=TRUE);   SPs=SPs[x$ix,];
    toRm=c();
    for  (sp in 1:nrow(SPs)){
      ii=sp; ##sp_i * x != sp_j for all x in 2:6 and all SP pairs (i,j)
      x=2; ##Check for doublets only
      maxDev=SPs[sp,'precision']* 2/3 ;
      i_=which(abs(SPs[sp,'Mean Weighted']*x-SPs[,'Mean Weighted'])<maxDev);
      ii=union(ii,i_); 
      if( length(ii)>1 ){
        ii=ii[which(SPs[ii,'nMutations']<0.6*max(SPs[ii,'nMutations'],na.rm=T))];; ##Keep only SP of max kurtosis
        toRm=c(toRm,ii )
      }
    }
    
    toRm=unique(toRm)
    if (!is.null(toRm) && length(toRm)>0){
      iReassign=which(dm[,'SP'] %in% SPs[toRm,'Mean Weighted'] | dm[,'SP_cnv'] %in% SPs[toRm,'Mean Weighted'] );
      print(paste('Reassigning SNVs after pruning',length(toRm),'doublet subpopulation(s). Pruned subpopulation(s):'))
      print(SPs[toRm,'Mean Weighted'])
      SPs=SPs[-toRm,, drop=FALSE];
      aM = assignMutations( dm[iReassign,,drop=FALSE], SPs, peakselection=peakselection);
      dm[iReassign,]=aM$dm; 
      iSP=match(aM$finalSPs[,'Mean Weighted'],SPs[,'Mean Weighted']);
      SPs[iSP,'nMutations']=SPs[iSP,'nMutations']+aM$finalSPs[,'nMutations'];
      finalSPs=SPs;
    }else{
      continueprune=FALSE;
      finalSPs=SPs;
    }
  }
  
  #if (plotF>2){
  #    if(!require(rgl)){
  #	message("Plot supressed: Package rgl required for 3D plot of subpopulation clusters. Load this package before using this option.")
  #    }else{
  #    ##plot probability distributions
  #    cols=c("red","yellow","green","pink","magenta","cyan","lightblue","blue");
  #    kk=ceil(nrow(finalSPs)/2); par(mfcol=c(2,kk));
  #    for (i in 2:nrow(finalSPs)){
  #        idx=which(dm[toUseIdx,"SP"]==finalSPs[i,"Mean Weighted"]);
  #        open3d();
  #	        persp3d(as.numeric(1:length(idx)),as.numeric(cfd$freq),t(cfd$densities[idx,]),col=cols[mod(i,length(cols))+1],aspect=c(1, 1, 0.5), add=FALSE,xlab="Mutation", ylab="cell-frequency", zlab="Probability");
  #	        title3d(paste("SP_", round(finalSPs[i,"Mean Weighted"], digits=2),"\noo",sep=""));
  #        	play3d(spin3d(axis=c(0,0,1), rpm=10), duration=5)
  #	    }
  #    }
  #}
  
  ##phylogeny
  finalSPs=.addColumn(finalSPs,'Ancestor',NA);
  finalSPs=.addColumn(finalSPs,'ClosestDescendant',NA);
  aQ=try(assignQuantityToSP(copyNumber, dm),silent=FALSE);
  tr=NULL;
  if(class(aQ)=="try-error" || is.null(ncol(aQ$ploidy))){
    print("Error encountered while assigning subpopulation specific ploidy to copy number segments. Phylogeny will not be inferred.")
  }else {
    dm=aQ$dm;
    aQ=aQ$ploidy;
    .writeExpandsOutput(X=aQ, dirF,snvF,suffix=".sps.cbs", message="Subpopulation specific copy numbers")
    output=paste(dirF, .Platform$file.sep, snvF, sep="");    # output=paste(dirF, .Platform$file.sep, gsub("\\.","_",snvF), sep="");
    tr=try(buildPhylo(aQ,output,dm=dm),silent=FALSE);
    if(class(tr)!="try-error" ){
      if(class(tr$dm)!="try-error" && !is.na(tr$dm)){
        dm=tr$dm;
        ##Assign ancestor and closest descendant
        phySPs=colnames(tr$spRelations)
        for (sp in rownames(tr$spRelations)){
          desc=as.numeric(gsub('SP_','',c(sp,phySPs[tr$spRelations[sp,]==1])))
          if(length(desc)>1){
            iAnc=which.min(abs(finalSPs[,'Mean Weighted']-desc[1]))
            desc=desc[which.min(abs(desc[2:length(desc)]-desc[1]))+1]
            iDesc=which.min(abs(finalSPs[,'Mean Weighted']-desc))
            finalSPs[iDesc,'Ancestor']=finalSPs[iAnc,'Mean Weighted']
            finalSPs[iAnc,'ClosestDescendant']=finalSPs[iDesc,'Mean Weighted']
          }
        }
      }else{
        print("Error encountered while reconstructing phylogeny")
      }
      tr=tr$tree;
    }
  }
  
  ####################
  ###Write output#####
  .writeExpandsOutput(X=dm, dirF,snvF,suffix=".sps", message="Subpopulation specific point mutations")
  .writeExpandsOutput(X=finalSPs, dirF,snvF,suffix=".spstats", message="Summary file of detected subpopulations")
  
  if (plotF>0){
    try(plotSPs(dm[toUseIdx,],snvF),silent=FALSE);
  }
  if (plotF>1 && !is.null(tr) && class(tr)!="try-error"){
    try(plot(tr),silent=FALSE);
  }
  
  result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities,"ploidy"=aQ,"tree"=tr,"homDelRegions"=homDelRegions);
  return(result);
}

.writeExpandsOutput<-function(X, dirF,snvF,suffix=".sps.cbs", message="Output"){
  output=paste(dirF, .Platform$file.sep, snvF,suffix,sep="");
  write(paste("## expands version",packageVersion("expands")), file = output, append=FALSE);
  suppressWarnings(write.table(X,file = output, append=TRUE, quote = FALSE, sep = "\t", row.names=FALSE));
  print(paste(message," saved under ",output));
}
