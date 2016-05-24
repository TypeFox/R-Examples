checkDNAbcd <-
function(seq,DistModel="K80") {   
 mysplit<-c();
 mylabels<-c();
 listsp<-data.frame();
 DNAlength<-c();
 mysplit<-strsplit(labels(seq),"_");
 mygen<-c(); for (i in 1:length(mysplit)){ mygen<-c(mygen,mysplit[[i]][1]);}
 mysp<-c(); for (i in 1:length(mysplit)){ mysp<-c(mysp,mysplit[[i]][2]);}
 mylabels<-data.frame(genus=mygen,species=mysp,id=labels(seq));
 listsp<-data.frame(species=sort(unique(paste(mylabels$genus,mylabels$species,sep="_"))),Nseq=NA,Nhap=NA);
 for(i in 1:length(listsp$species)) {
  listsp$Nseq[[i]]<-length(grep(listsp$species[[i]],labels(seq)));  
  if (listsp$Nseq[[i]]>1){
   listsp$Nhap[[i]]<-dim(haplotype(seq[grep(listsp$species[[i]],labels(seq)), ]))[[1]];  
  }
  else {
   listsp$Nhap[[i]]<- listsp$Nseq[[i]];
  }
 } 
 DNAlength<-dim(seq)[2]-as.numeric(checkDNA(seq));
 dist<-dist.dna(seq,DistModel,as.matrix=TRUE,pairwise.deletion=TRUE);
 diag(dist)<-NA
 spdist<-c();
 spdist<-sppDist(dist,paste(mylabels$genus,mylabels$species,sep="_"))
 colnames(dist)<-paste(mylabels$genus,mylabels$species,sep="_");
 write.csv(mylabels,"mylabels.csv"); 
 write.csv(listsp,"listsp.csv"); 
return(list(mylabels=mylabels, listsp=listsp, DNAlength=DNAlength, dist=dist, spdist=spdist, seq=seq));
}
