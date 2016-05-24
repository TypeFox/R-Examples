assignQuantityToMutation<-function(dm,cbs,quantityColumnLabel="CN_Estimate"){

##Add column segmentID to CBS
cols=c("quantityID",quantityColumnLabel,"segmentLength");
if (!any(colnames(cbs)==cols[1])){
    cbs=.addColumn(cbs,cols[1],NA);
    if (any(colnames(cbs)=="Count")){
        cbs[,"quantityID"]=cbs[,"Count"];
    }else{
        cbs[,"quantityID"]=t(1:nrow(cbs));
    }
}

##First add columns to input data if necessary
for (k in 1:length(cols)){
    dm=.addColumn(dm,cols[k],NA);
}
if (!any(colnames(cbs)=="segmentLength")){
    cbs=.addColumn(cbs,"segmentLength",NA);
    cbs[,"segmentLength"]=cbs[,"endpos"]-cbs[,"startpos"];
}

if (quantityColumnLabel[1]=="FPKM"){
    dm=.assignFPKMToMutation(dm,cbs,cols);
}else if(quantityColumnLabel[1]=="CN_Estimate" || !isempty(grep("^SP_",quantityColumnLabel[1])) ){
    dm=.assignCBSToMutation(dm,cbs,cols);
}else{
    stop(paste("Invalid quantityColumnLabel: ",quantityColumnLabel,". Valid options are: FPKM, CN_Estimate."));
}
return(dm);
}

.assignFPKMToMutation<-function(dm,cbs,cols){
dm=.addColumn(dm,"Dominant_Isoform",0);
dmPlus=c();
print("Assigning expression to mutations...")
##Assign copy numbers in cbs to mutations in dm
for (k in 1:nrow(dm)){
    if (mod(k,100)==0){
        print(paste("Finding overlaps for mutation", k, "out of ",nrow(dm),"..."))
    }
    idx=which(dm[k,"chr"]==cbs[,"chr"] & dm[k,"startpos"]>=cbs[,"startpos"] & dm[k,"startpos"]<=cbs[,"endpos"]);
    if (length(idx)>0){
        dm[k,cols]=cbs[idx[1],cols];
        dm[k,"Dominant_Isoform"]=dm[k,"Dominant_Isoform"]+1;
        idx=idx[-1]
    }
#    if (length(idx)>0){
#        dm2=repmat(dm[k,],length(idx),1);
#        dm2[,cols]=cbs[idx,cols];
#        dmPlus=rbind(dmPlus,dm2);
#    }
}
dm=rbind(dm,dmPlus);
dm=dm[,colnames(dm)!="segmentLength"];
print("... Done.")

return(dm);
}


.assignCBSToMutation<-function(dm,cbs,cols){
print("Assigning copy number to mutations...")
##Assign copy numbers in cbs to mutations in dm
for (k in 1:nrow(cbs)){
    if (mod(k,100)==0){
        print(paste("Finding overlaps for CBS segment", k,"out of ",nrow(cbs),"..."));
    }
    idx=which(dm[,"chr"]==cbs[k,"chr"] & dm[,"startpos"]>=cbs[k,"startpos"] & dm[,"startpos"]<=cbs[k,"endpos"]);
    if (length(idx)==0){
        next;
    }
    ok=which(is.na(dm[idx,"segmentLength"]) | dm[idx,"segmentLength"]>cbs[k,"segmentLength"]);
    if (length(ok)==0){
        next;
    }
    dm[idx[ok],cols]=repmat(cbs[k,cols],length(ok),1);
}
dm=dm[,colnames(dm)!="segmentLength"];
print("... Done.")

return(dm);
}
