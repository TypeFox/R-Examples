## Julia Bischof
## 10-09-2015

#library(doParallel)
#library(parallel)

clones<-function(aaseqtab=NULL,summarytab=NULL, ntseqtab=NULL,identity=0.85, useJ=TRUE,dispD=FALSE, 
                 dispSeqID=FALSE,dispCDR3aa=FALSE,dispCDR3nt=FALSE,dispJunctionFr.ratio=FALSE,
                 dispJunctionFr.list=FALSE,dispFunctionality.ratio=FALSE,dispFunctionality.list=FALSE,
                 dispTotalSeq=FALSE,nrCores=1){
  if(length(aaseqtab)==0 || length(grep("CDR3_IMGT",colnames(aaseqtab)))==0){
    stop("--> 5_AA-sequences file is missing")
  }
  if(length(summarytab)==0 && dispTotalSeq==T && length(grep("Sequence",colnames(summarytab)))==0){
    stop("--> 1_Summary file is missing, no information about Total_sequences_nt available")
  }
  if((length(summarytab)==0 && (dispJunctionFr.list==T || dispJunctionFr.ratio==T)) || (length(summarytab)>0 && length(grep("JUNCTION_frame",colnames(summarytab)))==0)){
    stop("--> 1_Summary file is missing, no information about JUNCTION_frames available")
  }
  if((length(ntseqtab)==0 && dispCDR3nt==T ) || (length(ntseqtab)>0 &&length(grep("CDR3_IMGT",colnames(ntseqtab)))==0)){
    stop("--> 3_Nt-sequences file is missing, no information about CDR3 nt sequences available")
  }
  if(as.numeric(nrCores)>as.numeric(detectCores())){
    stop(paste("--> nrCores is higher than available number of cores (only ",as.numeric(detectCores())," cores available)",sep=""))
  }
  
  cl<-makeCluster(nrCores)
  registerDoParallel(cl)
  
  V<-unlist(unique(apply(data.frame(aaseqtab$V_GENE_and_allele),1,function(x){strsplit(x,split=" |,|;|[*]")[[1]]})))
  V<-unique(V[grep("V",V)])
  J<-unlist(unique(apply(data.frame(aaseqtab$J_GENE_and_allele),1,function(x){strsplit(x,split=" |,|;|[*]")[[1]]})))
  J<-unique(J[grep("J",J)])
  
  tempout<-vector()
  i<-NULL
  clonelist<-foreach(i=1:length(V)) %dopar% {
    if(useJ==T){ ### useJ=T
      for(j in J){
        aaseqtab.sub<-aaseqtab[intersect(grep(paste(V[i],"[!/*]",sep=""),aaseqtab$V_GENE_and_allele,perl=T),grep(paste(j,"[!/*]",sep=""),aaseqtab$J_GENE_and_allele,perl=T)),c('CDR3_IMGT','V_GENE_and_allele','J_GENE_and_allele','D_GENE_and_allele','Functionality','Sequence_ID')]
        ntseqtab.sub<-data.frame(ntseqtab[intersect(grep(paste(V[i],"[!/*]",sep=""),aaseqtab$V_GENE_and_allele,perl=T),grep(paste(j,"[!/*]",sep=""),aaseqtab$J_GENE_and_allele,perl=T)),'CDR3_IMGT'])
        summarytab.sub<-summarytab[intersect(grep(paste(V[i],"[!/*]",sep=""),aaseqtab$V_GENE_and_allele,perl=T),grep(paste(j,"[!/*]",sep=""),aaseqtab$J_GENE_and_allele,perl=T)),c('JUNCTION_frame','Sequence')]
        
        CDR3length<-apply(data.frame(aaseqtab.sub$CDR3_IMGT),1,function(x){nchar(x)})
        uniqueCDR3length<-sort(as.numeric(unique(CDR3length)))
        if(length(uniqueCDR3length)==1){ ### only 1 CDR3 length
          temp<-vector()
          uniqueCDR3.sub<-unique(aaseqtab.sub$CDR3_IMGT)
          a.dist<-adist(unique(uniqueCDR3.sub))  
          a.dist[upper.tri(a.dist,diag = T)]<-NA
          tr<-floor(uniqueCDR3length*(1-identity))
          for(a in nrow(a.dist):1){
            if(length(which(a.dist[a,]<=tr))>0){
              if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
                temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))
                
                aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
                summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                
                tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
                                         uniqueCDR3length, # CDR3 length
                                         length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
                                         nrow(aaseqtab.new), # number all sequences, belonging to clone
                                         do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                         V[i], # V_gene
                                         do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                         if(useJ==T){j}else{do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", "))}, # J_gene
                                         do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                         
                                         if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                         if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                         if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                         if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                         if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                         if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                         if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                         if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                         if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
              }
            }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), aaseqtab.sub$CDR3_IMGT))>1){
              temp<-c(temp,uniqueCDR3.sub[a])
              
              aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
              ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
              summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
              
              tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
                                       l, # CDR3 length
                                       length(uniqueCDR3.sub[a]), # number_shared_CDR3
                                       nrow(aaseqtab.new), # number all sequences, belonging to clone
                                       do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                       V[i], # V_gene
                                       do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                       if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                       do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                       
                                       if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                       if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                       if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                       if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                       if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                       if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                       if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                       if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                       if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
            }
          }
        }else{ ### several CDR3 length
          for(l in uniqueCDR3length){
            temp<-vector()
            uniqueCDR3.sub<-unique(aaseqtab.sub$CDR3_IMGT[which(CDR3length==l)])
            a.dist<-adist(unique(uniqueCDR3.sub))  
            a.dist[upper.tri(a.dist,diag = T)]<-NA
            tr<-floor(l*(1-identity))
            for(a in nrow(a.dist):1){
              if(length(which(a.dist[a,]<=tr))>0){
                if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
                  temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))
                  
                  aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                  ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
                  summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                  
                  tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
                                           l, # CDR3 length
                                           length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
                                           nrow(aaseqtab.new), # number all sequences, belonging to clone
                                           do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                           V[i], # V_gene
                                           do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                           if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                           do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                           
                                           if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                           if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                           if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                           if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                           if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                           if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                           if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                           if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                           if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                           if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                           if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                           if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                           if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
                }
              }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), aaseqtab.sub$CDR3_IMGT))>1){
                temp<-c(temp,uniqueCDR3.sub[a])
                
                aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
                summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                
                tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
                                         l, # CDR3 length
                                         length(uniqueCDR3.sub[a]), # number_shared_CDR3
                                         nrow(aaseqtab.new), # number all sequences, belonging to clone
                                         do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                         V[i], # V_gene
                                         do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                         if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                         do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                         
                                         if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                         if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                         if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                         if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                         if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                         if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                         if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                         if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                         if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
              }
            }
          }
        }
      }
    }else{ ### useJ=F
      aaseqtab.sub<-aaseqtab[grep(paste(V[i],"[!/*]",sep=""),aaseqtab$V_GENE_and_allele,perl=T),c('CDR3_IMGT','V_GENE_and_allele','J_GENE_and_allele','D_GENE_and_allele','Functionality','Sequence_ID')]
      ntseqtab.sub<-data.frame(ntseqtab[grep(paste(V[i],"[!/*]",sep=""),aaseqtab$V_GENE_and_allele,perl=T),'CDR3_IMGT'])
      summarytab.sub<-summarytab[grep(paste(V[i],"[!/*]",sep=""),aaseqtab$V_GENE_and_allele,perl=T),c('JUNCTION_frame','Sequence')]
      
      CDR3length<-apply(data.frame(aaseqtab.sub$CDR3_IMGT),1,function(x){nchar(x)})
      uniqueCDR3length<-sort(as.numeric(unique(CDR3length)))
      if(length(uniqueCDR3length)==1){ ### only 1 CDR3 length
        temp<-vector()
        uniqueCDR3.sub<-unique(aaseqtab.sub$CDR3_IMGT)
        a.dist<-adist(unique(uniqueCDR3.sub))  
        a.dist[upper.tri(a.dist,diag = T)]<-NA
        tr<-floor(uniqueCDR3length*(1-identity))
        for(a in nrow(a.dist):1){
          if(length(which(a.dist[a,]<=tr))>0){
            if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
              temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))
              
              aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
              ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
              summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
              
              tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
                                       uniqueCDR3length, # CDR3 length
                                       length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
                                       nrow(aaseqtab.new), # number all sequences, belonging to clone
                                       do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                       V[i], # V_gene
                                       do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                       if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                       do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                       
                                       if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                       if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                       if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                       if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                       if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                       if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                       if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                       if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                       if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
            }
          }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), aaseqtab.sub$CDR3_IMGT))>1){
            temp<-c(temp,uniqueCDR3.sub[a])
            
            aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
            ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
            summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
            
            tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
                                     l, # CDR3 length
                                     length(uniqueCDR3.sub[a]), # number_shared_CDR3
                                     nrow(aaseqtab.new), # number all sequences, belonging to clone
                                     do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                     V[i], # V_gene
                                     do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                     if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                     do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                     
                                     if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                     if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                     if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                     if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                     if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                     if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                     if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                     if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                     if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                     if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                     if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                     if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                     if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
          }
        }
      }else{ ### several CDR3 length
        for(l in uniqueCDR3length){
          temp<-vector()
          uniqueCDR3.sub<-unique(aaseqtab.sub$CDR3_IMGT[which(CDR3length==l)])
          a.dist<-adist(unique(uniqueCDR3.sub))  
          a.dist[upper.tri(a.dist,diag = T)]<-NA
          tr<-floor(l*(1-identity))
          for(a in nrow(a.dist):1){
            if(length(which(a.dist[a,]<=tr))>0){
              if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
                temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))
                
                aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
                summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
                
                tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
                                         l, # CDR3 length
                                         length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
                                         nrow(aaseqtab.new), # number all sequences, belonging to clone
                                         do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                         V[i], # V_gene
                                         do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                         if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                         do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                         
                                         if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                         if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                         if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                         if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                         if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                         if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                         if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                         if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                         if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                         if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
              }
            }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), aaseqtab.sub$CDR3_IMGT))>1){
              temp<-c(temp,uniqueCDR3.sub[a])
              
              aaseqtab.new<-aaseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
              ntseqtab.new<-data.frame(ntseqtab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),])
              summarytab.new<-summarytab.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",aaseqtab.sub$CDR3_IMGT),perl=T),]
              
              tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
                                       l, # CDR3 length
                                       length(uniqueCDR3.sub[a]), # number_shared_CDR3
                                       nrow(aaseqtab.new), # number all sequences, belonging to clone
                                       do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",aaseqtab.new$CDR3_IMGT),perl=T))})), sep=", ")), # sequence count
                                       V[i], # V_gene
                                       do.call(paste, c(as.list(unique(aaseqtab.new$V_GENE_and_allele)), sep=", ")), # V_gene & allele
                                       if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                                       do.call(paste, c(as.list(unique(aaseqtab.new$J_GENE_and_allele)), sep=", ")), # J_gene & allele
                                       
                                       if(dispD==T){do.call(paste, c(as.list(unique(aaseqtab.new$D_GENE_and_allele)), sep=", "))},
                                       if(dispCDR3aa==T){do.call(paste, c(as.list(aaseqtab.new$CDR3_IMGT), sep=", "))},
                                       if(length(ntseqtab.new)>0 && dispCDR3nt==T){do.call(paste, c(as.list(ntseqtab.new[,1]), sep=", "))},
                                       if(dispFunctionality.list==T){do.call(paste, c(as.list(aaseqtab.new$Functionality), sep=", "))},
                                       if(dispFunctionality.ratio==T){length(grep("^productive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                       if(dispFunctionality.ratio==T){length(grep("^unproductive",aaseqtab.new$Functionality))/nrow(aaseqtab.new)},
                                       if(dispFunctionality.ratio==T){length(grep("prod",aaseqtab.new$Functionality,invert=T))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.list==T){do.call(paste, c(as.list(summarytab.new$JUNCTION_frame), sep=", "))},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="in-frame"))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(which(summarytab.new$JUNCTION_frame=="out-of-frame"))/nrow(aaseqtab.new)},
                                       if(length(summarytab.new)>0 && dispJunctionFr.ratio==T){length(grep("frame",summarytab.new$JUNCTION_frame,invert=T))/nrow(aaseqtab.new)},                             
                                       if(dispSeqID==T){do.call(paste, c(as.list(aaseqtab.new$Sequence_ID), sep=", "))},
                                       if(length(summarytab.new)>0 && dispTotalSeq==T){do.call(paste, c(as.list(unique(summarytab.new$Sequence)), sep=", "))}))
            }
          }
        }
      }
    }
    return(data.frame(tempout,stringsAsFactors = F))
  }
  
  
  clonRel<-do.call(rbind.data.frame, clonelist)
  
  stopCluster(cl)
  
  if(is.data.frame(clonRel) && nrow(clonRel)>0){
    colnames(clonRel)<-c("unique_CDR3_sequences_AA",
                         "CDR3_length_AA",
                         "number_of_unique_sequences",
                         "total_number_of_sequences",
                         "sequence_count_per_CDR3",
                         "V_gene",
                         "V_gene_and_allele",
                         "J_gene",
                         "J_gene_and_allele",
                         if(dispD==T){"D_gene"},
                         if(dispCDR3aa==T){"all_CDR3_sequences_AA"},
                         if(length(ntseqtab)>0 && dispCDR3nt==T){"all_CDR3_sequences_nt"},
                         if(dispFunctionality.list==T){"Functionality_all_sequences"},
                         if(dispFunctionality.ratio==T){c("Func_productive_sequences","Func_unproductive_sequences", "Func_unknown")},
                         if(length(summarytab)>0 && dispJunctionFr.list==T){"Junction_frame_all_sequences"}, 
                         if(length(summarytab)>0 && dispJunctionFr.ratio==T){c("JF_in_frame","JF_out_of_frame","JF_unknown")},
                         if(dispSeqID==T){"Sequence_IDs"},
                         if(length(summarytab)>0 && dispTotalSeq==T){"Total_sequences_nt"})    
    
    clonRel<-clonRel[duplicated(clonRel[,"unique_CDR3_sequences_AA"])==F,]
    if(length(intersect(which(clonRel[,"total_number_of_sequences"]==1),
                        which(clonRel[,"sequence_count_per_CDR3"]==1)))>0){
      clonRel<-clonRel[-intersect(which(clonRel[,"total_number_of_sequences"]==1),
                                  which(clonRel[,"sequence_count_per_CDR3"]==1)),]
    }
  }else{
    clonRel<-NULL
    cat("... no clones\n")
  } 
  return(data.frame(clonRel,check.names=F,stringsAsFactors = F))  
}

