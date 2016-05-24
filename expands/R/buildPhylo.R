buildPhylo<-function(ploidy,outF,treeAlgorithm="bionjs",dm=NA, add="Germline"){
  #library(ape)
  out=list("tree"=NULL,"dm"=dm);
  ii=grep("SP",colnames(ploidy));
  cnv=ploidy[,ii];
  if(is.null(ncol(cnv))){
    print("Less than two SPs coexist in this tumor. Aborting phylogeny reconstruction");
    return(out);
  }
  
  print(paste("Building phylogeny using ",treeAlgorithm," algorithm",sep=""))
  print("Pairwise SP distances calculated as: % segments with identical copy number");
  ##Add user specified artificial SP, if any
  if(!is.null(add)){
    if (!any(colnames(cnv)==add,na.rm=T) && (add=="Consensus" || add=="Germline") ){
      cnv=cbind(cnv,matrix(matrix(NaN,nrow(cnv),1),nrow=nrow(cnv), ncol=1, dimnames=list(rownames(cnv),add)));
      if(add=="Consensus"){
        cnv[,add]=round(colMeans(t(cnv),na.rm=T)); ####ploidy:=consensus at all positions
      }else if(add=="Germline"){
        cnv[,add]=2; ##assuming diploid germline status
      }
    }
  }
  
  toRm=c(); 
  ##distance matrix from pairwise alignments
  cols=gsub(" ","",colnames(cnv));
  D=matrix(matrix(1,ncol(cnv),ncol(cnv)), nrow=ncol(cnv), ncol=ncol(cnv), dimnames=list(cols,cols));
  for (i in 1:ncol(cnv)){
    for (j in i:ncol(cnv)){
      ii=which(!is.na(cnv[,i]) & !is.na(cnv[,j]));
      if (length(ii)==0){
        D[i,j]<-D[j,i]<-NA;
        next;
      }
      x=cnv[ii,i];      y=cnv[ii,j];
      dd=length(which(x!=y))/length(ii);
      if (i!=j){
        dd=dd+0.3;
      }
      D[i,j]<-D[j,i]<-dd
    }
    if (any(is.na(D[i,1:i]))){
      toRm=cbind(toRm,i);#remove NAs
    }
  }
  #remove NAs
  if (length(toRm)>0){
    print(paste("Insufficient copy number segments for ",rownames(D)[toRm],". SP excluded from phylogeny",sep=""))
    D=D[-toRm,];
    D=D[,-toRm];
  }
  
  if(is.null(nrow(D)) || nrow(D)<3){
    print("No two SPs found between which distance could be calculated. Aborting phylogeny reconstruction");
    return(out);
  }
  
  D=D*100;
  write.table(D,paste(outF,".dist",sep=""),quote=F,sep="\t");
  print(paste("distance-matrix saved under ",outF,".dist",sep=""));
  tr =c();
  if (treeAlgorithm=="bionjs"){
    tr <- bionjs(D);
  }else{
    tr <- njs(D);
  }
  tr$root.edge <- 0; ## adds root dummy
  
  outF=paste(outF,".tree",sep="");
  write.tree(tr, file = outF);
  print(paste("tree saved under ",outF,sep=""));
  
  out$tree=tr;
  out$spRelations=NULL;
  if(!is.na(dm)){
    out1=try(.assignSNVsToMultipleSPs(dm,outF),silent=FALSE)
    if(class(out1)!="try-error"){
      dm=out1$dm;
      out$spRelations=out1$spRelations;
    }
    out$dm=dm;
  }
  return(out);
}

.assignSNVsToMultipleSPs <-function(dm,outF){
  if (!requireNamespace("phylobase", quietly = TRUE)) {
    print("Package \'phylobase\' is needed for assigning SNVs to Multiple SPs but is not installed.")
    return(dm);
  }
  tr=try(phylobase::readNewick(outF,check.names=F),silent=F);
  if(class(tr)=="try-error"){
    print("Warning! Setting negative edge lengths to 0!")
    tr=read.tree(outF);
    tr$edge.length[tr$edge.length<0]=0
    write.tree(tr, file = "tmp.tree");
    tr=phylobase::readNewick("tmp.tree",check.names=F);
  }
  print("Assigning SNVs to SPs...")
  spsInTree=names(phylobase::getNode(tr,type="tip"));  
  SPs = sort(unique(c(dm[, "SP_cnv"],dm[,"SP"])))
  spSizes = unique(round(SPs * 1000)/1000)
  spNames= paste("SP_", as.character(spSizes), sep = "");
  x = colnames(dm)
  x[(length(x) + 1):(length(x) + length(SPs))] =spNames
  x=c(x,"Clone");
  dm = cbind(dm, matrix(0, nrow(dm), length(SPs)), dm[,"SP"]);
  colnames(dm)=x;
  ##Save ancestor-to-descendant (rows-to-columns) relations:
  spRelations=matrix(0,length(spNames),length(spNames));
  rownames(spRelations)=spNames; colnames(spRelations)=spNames;
    
  for (k in 1:nrow(dm)) {
    thisSP = paste("SP_", as.character(round(dm[k,"SP"] * 1000)/1000),sep="");
    if (!is.na(dm[k,"SP"])){
      dm[k,gsub(" ","_",thisSP)]=1; #dm[k,"PM_B"]; binary assignment for now 
    }
    if(!thisSP %in% spsInTree){
      next;
    }
    if (mod(k, 100) == 0) {
      print(paste("Assigning SPs for SNV", k, 
                  "out of ", nrow(dm), "..."))
    }
    dm=.propagateSNVToMultipleSPs(thisSP,dm,k,tr,spSizes)
    spRelations[thisSP,dm[k,colnames(spRelations)]==1]=1;
    spRelations[thisSP,thisSP]=0
  }
  
  for (i in 1:(length(spNames)-1)){
    iPhylo=which(sum(t(dm[,spNames]!=0))>length(spNames)-i)
    print(paste(length(iPhylo), " SNVs assigned to >",length(spNames)-i," SPs"))
  }
  return(list(dm=dm,spRelations=spRelations) )
}





.propagateSNVToMultipleSPs <-function(thisSP,dm,k,tr, spSizes){
  thisSPsize=as.numeric(gsub("SP_","",thisSP));
  xx=phylobase::getNode(tr,type="tip");
  ii_This=match(thisSP,names(xx)); ##Node representing SP which harbors this SNV
  sibl=phylobase::siblings(tr,xx[ii_This]); ##Siblings of SP with this SNV
  if (length(sibl)==0){
    return(dm)
  }
  for (s in 1:length(sibl)){
    desc=.getAllDescendingTips_InclSelf(tr,sibl[s],c());
    i_toKEEP=union(grep("SP",names(desc)),which(is.na(names(desc))));
    desc=desc[i_toKEEP];
    for (sd in names(desc)){
      otherSP=gsub(" ","_",sd);
      otherSPsize=as.numeric(gsub("SP_","",otherSP))
      rootL=0
      if(is.na(names(sibl[s])) || sd!=names(sibl[s])){
        rootL=phylobase::edgeLength(tr,sibl[s]);
      }
      if (dm[k,otherSP]==0 && rootL+ phylobase::edgeLength(tr,sd) > 1.25*phylobase::edgeLength(tr,xx[ii_This])){ ##1.4
        ##This SP is likely ancestor of SP assigned as sibling. TODO: find tree reconstruction algorithm that can assign "living populations" as common ancestors
        dm[k,otherSP]=1; #dm[k,"PM_B"];  ##Other SP inherits mutation of this SP. 
        dm[k,"Clone"]=thisSPsize-otherSPsize;
        ##TODO: what if the ploidy of the otherSP in this region is < B-allele ploidy of this SP! need to check this and set to minimum. Temporary solution --> binary assignment (above)
        dm=.propagateSNVToMultipleSPs(otherSP,dm,k,tr,spSizes)
      }
    }
  }
  return(dm);
}




.getAllDescendingTips_InclSelf<-function(tr,nd,desc){
  kids=phylobase::children(tr, nd);
  if(length(kids)==0){
    return(nd)
  }
  for (i in 1:length(kids)){
    if(is.na(names(kids)[i])){
      desc=c(desc,.getAllDescendingTips_InclSelf(tr,kids[i],desc));
    }else{
      desc=c(desc,kids[i]);
    }
  }
  return(desc);
}

