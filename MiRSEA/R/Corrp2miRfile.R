 ###############################################
 ##get p-value from miR and pathway , p2miR(pathway-miR file)
 Corrp2miRfile<-function(pathway="kegg",species="hsa"){
   pathway <- GetPathwayData(pathway)
   MiRTarget <- GetMiRTargetData()
   if(species=="hsa"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="hsa"),]
   }else if(species=="ath"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="ath"),]
   }else if(species=="cel"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="cel"),]
   }else if(species=="dre"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="dre"),]
   }else if(species=="dme"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="dme"),]
   }else if(species=="mmu"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="mmu"),]
   }else if(species=="rno"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="rno"),]
   }else if(species=="example"){
   file_spe<-MiRTarget[which(MiRTarget[,4]=="hsa"),]
   file_spe<-file_spe[1:200,]
   }
   
   spe.file<-file_spe[,2:3];
   spe.file<-unique(spe.file);
   m.file<-spe.file[,1]
   m.file <-unique(m.file);
   mm <- vector(length = length(m.file), mode = "numeric") 
   nn <- vector(length = length(pathway), mode = "numeric")
   max.m <- 0	
   pathway.names <- vector(length = length(pathway), mode = "character")
   pathway.dec <- vector(length = length(pathway), mode = "character")
   miR.names <- vector(length = length(m.file), mode = "character")
   p <- matrix(0, nrow=length(pathway), ncol= length(m.file))
     

     for (i in 1:length(pathway)) {
        #print(i)
        ps.line <- noquote(unlist(strsplit(pathway[[i]], "\t")));
        ps.line<-ps.line[ps.line!=""]
        pathway.names[i]<-ps.line[1]
        pathway.dec[i]<-ps.line[2]
        ps.line<-ps.line[-(1:2)]
        nn[i]<-length(ps.line)
        for(j in 1:length(m.file)){
             z <- spe.file[which(spe.file[,1]==m.file[j]),]
             miR <- vector(length = nrow(z), mode = "character")
              miR.names[j]<-as.character(z[1,1])
             for(r in 1:nrow(z)){
                 miR[r]<-as.character(z[r,2])
                }

            mm[j]<-length(miR)
        k <- length(which(miR%in%ps.line==TRUE))
		p[i,j]<-1-phyper(k-1,mm[j],26232-mm[j],nn[i])
         }
		max.m<-max(max.m,length(which(as.numeric(p[i,])<1))) 
     }
    p2m.count <- 1
    p2miR <- matrix("",nrow=length(pathway),ncol=max.m+2)
	 for(i in 1:nrow(p)){
        if(length(which(as.numeric(p[i,])<1))>0){
             p.miR<-miR.names[which(as.numeric(p[i,])<1)]
             p.miR<-append(pathway.dec[i],p.miR)
             p.miR<-append(pathway.names[i],p.miR)
             p2miR[p2m.count,1:length(p.miR)] <- p.miR 
	         p2m.count<- p2m.count+1
	    }
    }
	rownames(p)<-pathway.names
	colnames(p)<-miR.names
 return(list(p = p,p2miR = p2miR))
}