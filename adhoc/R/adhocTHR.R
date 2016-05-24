adhocTHR <-
function(a,NbrTh=30,ErrProb=0.05,Ambig="incorrect",Reg="linear") {   
 myBM<-list();#will gather the labels of the BM (species names)
 myBMid<-list();#will gather the indices of the BM
 myth<-c();

#### BM ANALYSIS ####
BM<-data.frame(labels=a$mylabels$id,distBM=NA,idBM=NA,IDcheck=NA);
 for(i in 1:length(a$mylabels$id)){
  spname<-c();
  BM$distBM[i]<-min(a$dist[i,],na.rm=TRUE);
  myBM[[i]]<-labels(which(a$dist[i,]==min(a$dist[i,],na.rm=TRUE)));   
  myBMid[[i]]<-labels(which(a$dist[,i]==min(a$dist[,i],na.rm=TRUE))); 
  spname<-paste(a$mylabels$genus[i],a$mylabels$species[i],sep="_");               
  if (length(unique(c(spname,myBM[[i]])))>1) {#When more than one species name is found among the best matches...
   BM$idBM[i]<-unique(c(spname,myBM[[i]]))[-which(unique(c(spname,myBM[[i]]))==spname)][[1]];  #...only one heterospecific allospecific name is given in idBM
  }
  else {
   BM$idBM[i]<-myBM[[i]][1];
  }
 }
 
#### BCM ANALYSIS ####
myth<-seq(from=0, to=as.numeric(max(BM$distBM,na.rm=TRUE)),length=NbrTh);
 allmatches<-data.frame(labels=a$mylabels$id);
 for(i in 1:length(a$mylabels$id)){
  for(j in 1:length(myth)){
   if (BM$distBM[i]>myth[j]){  #BM above threshold
    if (length(intersect(paste(a$mylabels$genus[i],a$mylabels$species[i],sep="_"),myBM[[i]]))==0) {
     allmatches[i,j+1]<-"TN";
    }
    else {
     if (length(unique(c(paste(a$mylabels$genus[i],a$mylabels$species[i],sep="_"),myBM[[i]])))>1) {
      allmatches[i,j+1]<-"TNambiguous" 
     }
     else{
      allmatches[i,j+1]<-"FN";
     }
    }
   }   
   else {  #BM below threshold
    if (length(intersect(paste(a$mylabels$genus[i],a$mylabels$species[i],sep="_"),myBM[[i]]))==0) {
     allmatches[i,j+1]<-"FP";
    }
    else {
     if (length(unique(c(paste(a$mylabels$genus[i],a$mylabels$species[i],sep="_"),myBM[[i]])))>1) {
      allmatches[i,j+1]<-"FPambiguous" 
     }
     else{
      allmatches[i,j+1]<-"TP";
     }
    }
   }
  }
 }
 BM$IDcheck<-allmatches[,length(myth)+1];
 colnames(allmatches)[2:(length(myth)+1)]<-paste("threshold nb",1:length(myth));

#### SUMMARY AND CALCULATION OF RE, OE, ETC ####
 IDcheck<-data.frame(thres=myth, TP=NA, FP=NA, RE=NA, TN=NA, FN=NA, OE=NA, Accuracy=NA, Precision=NA);
 for (k in 1:length(myth)) {
  if (Ambig=="correct"){
   IDcheck$TP[k]<-length(which(allmatches[,k+1]=="TP"))+length(which(allmatches[,k+1]=="FPambiguous"));
   IDcheck$FP[k]<-length(which(allmatches[,k+1]=="FP"));
   IDcheck$TN[k]<-length(which(allmatches[,k+1]=="TN"));
   IDcheck$FN[k]<-length(which(allmatches[,k+1]=="FN"))+length(which(allmatches[,k+1]=="TNambiguous"));
  }
  else {
   if (Ambig=="incorrect"){
    IDcheck$TP[k]<-length(which(allmatches[,k+1]=="TP"));
    IDcheck$FP[k]<-length(which(allmatches[,k+1]=="FP"))+length(which(allmatches[,k+1]=="FPambiguous"));
    IDcheck$TN[k]<-length(which(allmatches[,k+1]=="TN"))+length(which(allmatches[,k+1]=="TNambiguous"));
    IDcheck$FN[k]<-length(which(allmatches[,k+1]=="FN"));
   }
   else {
    IDcheck$TP[k]<-length(which(allmatches[,k+1]=="TP"));
    IDcheck$FP[k]<-length(which(allmatches[,k+1]=="FP"));
    IDcheck$TN[k]<-length(which(allmatches[,k+1]=="TN"));
    IDcheck$FN[k]<-length(which(allmatches[,k+1]=="FN"));
   }
  } 
  IDcheck$RE[k]<-IDcheck$FP[k]/(IDcheck$TP[k]+IDcheck$FP[k]);
  IDcheck$OE[k]<-(IDcheck$FP[k]+IDcheck$FN[k])/(IDcheck$TP[k]+IDcheck$FP[k]+IDcheck$TN[k]+IDcheck$FN[k]);
  IDcheck$Accuracy[k]<-(IDcheck$TP[k]+IDcheck$TN[k])/(IDcheck$TP[k]+IDcheck$FP[k]+IDcheck$TN[k]+IDcheck$FN[k]);
  IDcheck$Precision[k]<-IDcheck$TP[k]/(IDcheck$TP[k]+IDcheck$FP[k]);
 }
 write.csv(BM,"Bestmatches.csv"); 
 write.csv(IDcheck,"ID.csv"); 

#### REDFLAGGED SPECIES ####
 redflagged<-list();
 redflaggedSP<-c();
 myBMid<-lapply(myBMid,paste,collapse=" ");
 if (length(which(BM$IDcheck=="FPambiguous"))>0){#if there are ambiguous ID
  for (i in 1:length(myBMid[which(BM$IDcheck=="FPambiguous")])){
   redflagged[[i]]<-as.character(BM$labels[which(BM$IDcheck=="FPambiguous")][[i]]);#ID of the flagged sequence
   redflagged[[i]][2]<-length(unique(as.character(myBM[which(BM$IDcheck=="FPambiguous")][[i]])));#number of species in the complex
   redflagged[[i]][3]<-length(which(as.character(myBM[which(BM$IDcheck=="FPambiguous")][[i]])==paste(a$mylabels$genus[which(BM$IDcheck=="FPambiguous")][[i]],a$mylabels$species[which(BM$IDcheck=="FPambiguous")][[i]],sep="_")));#number of conspecific sequences
   redflagged[[i]][4]<-length(myBM[which(BM$IDcheck=="FPambiguous")][[i]])-as.numeric(redflagged[[i]][3]);#number of allospecific sequences
   redflagged[[i]][5]<-as.character(myBMid[which(BM$IDcheck=="FPambiguous")][[i]]);#listing of all Best matches
   redflaggedSP[[i]]<-unique(as.character(myBM[which(BM$IDcheck=="FPambiguous")][[i]]));
  }
  writeLines(c("red_flagged_seqID,Nb_species,Nb_conspecific_seq,Nb_allospecific_seq,list_of_best_matches",unlist(lapply(redflagged,paste,collapse=","))),"redflagged.csv"); 
  redflaggedSP<-unique(redflaggedSP);
  writeLines(c("red_flagged_species_groups",unlist(lapply(redflaggedSP,paste,collapse=","))),"redflaggedSP.csv"); 
 }

#### LINEAR REGRESSION ####
  if (Reg=="linear") {
   myreg<-c();
   THR<-c();
   myreg<-lm(RE~thres,IDcheck);
   myreg$coefficients[3:4]<-0;# this is to uniformise the format of the regression with the polynomial regression (coefficient for x^2 and x^3 = 0)
  }

#### POLYNOMIAL REGRESSION ####
  if (Reg=="polynomial") {
   myreg<-c();
   fp<-c();
   solp<-c();
   THRp<-c();
   IDcheck$thres2<-IDcheck$thres^2;
   IDcheck$thres3<-IDcheck$thres^3;
   myreg<- lm(RE ~ thres + thres2 + thres3, IDcheck);
   fp<-polynomial(myreg$coefficient);
   solp<-solve(fp,ErrProb);
   THR<-solp[(solp > 0) & (solp < max(IDcheck$thres))];
  }
   
#### OUTPUT ####
  if (length(grep("FP", BM$IDcheck))==0) stop("All identifications are correct when using the best match method (no distance threshold considered). An ad hoc distance threshold for best close match identification cannot be calculated");
  THR<-(ErrProb-as.numeric(myreg$coefficient[1]))/as.numeric(myreg$coefficient[2]);
  if (THR<0) stop("The estimated relative identification error (RE) cannot be reached using this reference library");
  return(list(BM=BM,IDcheck=IDcheck, reg=myreg, ErrProb=ErrProb, THR=THR, redflagged=redflagged, redflaggedSP=redflaggedSP));
}
