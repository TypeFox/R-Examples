## Julia Bischof
## 11-24-2015

#library(doParallel)
#library(parallel)


clones.shared<-function(clones.tab=NULL,identity=0.85,useJ=TRUE,dispD=TRUE,dispCDR3aa=FALSE,dispCDR3nt=FALSE,dispFunctionality.list=FALSE,
                        dispFunctionality.ratio=FALSE,
                        dispJunctionFr.list=FALSE,dispJunctionFr.ratio=FALSE,dispTotalSeq=FALSE,nrCores=1){
  if(length(clones.tab)==0){
    stop("--> Clone file is missing or empty")
  }
  if(as.numeric(nrCores)>as.numeric(detectCores())){
    stop(paste("--> nrCores is higher than available number of cores (only ",as.numeric(detectCores())," cores available)",sep=""))
  }
  
  uniCDR3length<-sort(unique(clones.tab$CDR3_length_AA))
    
  cl<-makeCluster(nrCores)
  registerDoParallel(cl)
  
  clones.uniV<-unique(unlist(apply(data.frame(clones.tab$V_gene),1,function(x){strsplit(x,split=" |,|;|_")[[1]]})))
  clones.uniV<-clones.uniV[grep(substr(clones.tab$V_gene[1],1,4),clones.uniV)]
  clones.uniV.gene<-unique(apply(data.frame(clones.uniV),1,function(x){strsplit(x,split="[*]")[[1]][1]}))
  
  if(useJ==TRUE){
    clones.J<-apply(data.frame(clones.tab$J_gene),1,function(x){strsplit(x,split=" |,|;|_")[[1]]})
    clones.uniJ<-unique(unlist(apply(data.frame(unlist(clones.J)),1,function(x){strsplit(x,split=" ")[[1]]})))
    clones.uniJ<-clones.uniJ[grep(substr(clones.tab$J_gene[1],1,4),clones.uniJ)]
    clones.uniJ.gene<-unique(apply(data.frame(clones.uniJ),1,function(x){strsplit(x,split="[*]")[[1]][1]}))
  }
  
  sharedclone.temp<-vector()
  i<-NULL
  clonelist<-foreach(i=1:length(uniCDR3length)) %dopar%{
    #for(i in 1:length(uniCDR3length)){
    #print(paste("... Analyzing clone with CDR3 length = ",uniCDR3length[i]," AA (",i,"/",length(uniCDR3length),")",sep=""))
    CDR3.index<-which(clones.tab$CDR3_length_AA==uniCDR3length[i])
    if(length(unique(CDR3.index))>1){
      for(v in 1:length(clones.uniV.gene)){
        V.index<-which(clones.tab$V_gene==clones.uniV.gene[v])                
        if(useJ==TRUE){ # use J
           for(j in 1:length(clones.uniJ.gene)){
             J.index<-which(clones.tab$J_gene==clones.uniJ.gene[j])
             CDR3VJ.index<-intersect(intersect(V.index,J.index), CDR3.index)
            if(length(CDR3VJ.index)>0){
              #cat(clones.uniV.gene[v],", ", clones.uniJ.gene[j],"... ")
              tr<-floor(uniCDR3length[i]*(1-as.numeric(identity)))
              CDR3.adist.temp<-unlist(lapply(clones.tab$unique_CDR3_sequences_AA[CDR3VJ.index],function(x){strsplit(x,split=" |,|;|_|-|[|]")[[1]]}))
              CDR3.adist<-unique(unlist(lapply(clones.tab$unique_CDR3_sequences_AA[CDR3VJ.index],function(x){strsplit(x,split=" |,|;|_|-|[|]")[[1]]})))
              CDR3.adist<-CDR3.adist[which(CDR3.adist!="")]
              ind.adist<-unlist(apply(data.frame(clones.tab[CDR3VJ.index,]),1,function(x){rep(x[1],length(strsplit(x[2],split=" |,|;|_|-|[|]")[[1]]))}))
              for(a in 1:length(CDR3.adist)){
                a.dist<-as.logical(adist(CDR3.adist[a], CDR3.adist)<=tr)
                a.dist[a]<-NA
                clones.tab.new<-vector()
                grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(unlist(apply(expand.grid(c(" ","^"),as.character(unique(CDR3.adist[c(a,which(a.dist==T))])), c(",","$")),1,paste, collapse=""))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
                if(length(which(a.dist==T))>0 && length(unique(ind.adist[grepCDR3]))>1){
                  #print(a)
                  clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(unlist(apply(expand.grid(c(" ","^"),as.character(unique(CDR3.adist[c(a,which(a.dist==T))])), c(",","$")),1,paste, collapse=""))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                  nrseq.temp<-vector()
                  cdr3temp<-unique(c(unlist(apply(data.frame(clones.tab.new$unique_CDR3_sequences_AA),1,function(x){strsplit(x, split=", ")[[1]]}))))
                  for(c in 1:length(cdr3temp)){
                    cdr3tempnr<-vector()
                    for(g in 1:length(unique(ind.adist[grepCDR3]))){
                      inter<-intersect(grep(unique(ind.adist[grepCDR3])[g],clones.tab.new$samples),grep(gsub("[*]","_",cdr3temp[c]),gsub("[*]","_",clones.tab.new$all_CDR3_sequences_AA), perl=T))
                      if(length(inter)==0){
                        cdr3tempnr<-c(cdr3tempnr, paste(unique(ind.adist[grepCDR3])[g], 0, sep=": "))
                      }else{
                        cdr3tempnr<-c(cdr3tempnr, paste(unique(ind.adist[grepCDR3])[g], length(gregexpr(gsub("[*]","_",cdr3temp[c]),gsub("[*]","_",c(clones.tab.new$all_CDR3_sequences_AA[inter])))[[1]])/length(inter), sep=": "))
                      }
                    }
                    nrseq.temp<-c(nrseq.temp, paste(cdr3temp[g]," (",do.call(paste, c(as.list(as.character(cdr3tempnr)), sep=", ")),")",sep=""))                   
                  }
                  nrseq<-do.call(paste, c(as.list(as.character(nrseq.temp)), sep="; "))
                  if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                    #print(paste(clones.uniV.gene[v], clones.uniJ.gene[j],a))
                    sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                               do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                               uniCDR3length[i],                                                               
                                                               do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist==T))]))), sep="; ")),
                                                               length(unique(CDR3.adist[c(a,which(a.dist==T))])),
                                                               nrseq,
                                                               do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                               clones.uniV.gene[v],
                                                               if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                               
                                                               if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                               if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                               if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                               if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                               if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                               if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                  } ## end: new entry
                  rm(clones.tab.new)
                } ## end if shared CDR3 exists
                rm(a.dist)
              } ## end for adist
            } ## end CDR3VJ.index>0
            rm(CDR3VJ.index)
            rm(CDR3.adist.temp)
            rm(CDR3.adist)
            rm(ind.adist)
            gc()
          } ## end J genes
        }else{ # use no J
          tr<-floor(uniCDR3length[i]*(1-as.numeric(identity)))
          
          CDR3V.index<-intersect(CDR3.index, V.index)
          CDR3.adist.temp<-unlist(lapply(clones.tab$unique_CDR3_sequences_AA[CDR3V.index],function(x){strsplit(x,split=" |,|;|_|-|[|]")[[1]]}))
          CDR3.adist<-unique(unlist(lapply(clones.tab$unique_CDR3_sequences_AA[CDR3V.index],function(x){strsplit(x,split=" |,|;|_|-|[|]")[[1]]})))
          CDR3.adist<-CDR3.adist[which(CDR3.adist!="")]
          ind.adist<-unlist(apply(data.frame(clones.tab[CDR3V.index,]),1,function(x){rep(x[1],length(strsplit(x[2],split=" |,|;|_|-|[|]")[[1]]))}))
          for(a in 1:length(CDR3.adist)){
            a.dist<-as.logical(adist(CDR3.adist[a], CDR3.adist)<=tr)
            a.dist[a]<-NA
            clones.tab.new<-vector()
            grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(unlist(apply(expand.grid(c(" ","^"),as.character(unique(CDR3.adist[c(a,which(a.dist==T))])), c(",","$")),1,paste, collapse=""))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
            if(length(which(a.dist==T))>0 && length(unique(ind.adist[grepCDR3]))>1){
              clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(unlist(apply(expand.grid(c(" ","^"),as.character(unique(CDR3.adist[c(a,which(a.dist==T))])), c(",","$")),1,paste, collapse=""))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
              nrseq.temp<-vector()
              cdr3temp<-unique(c(unlist(apply(data.frame(clones.tab.new$unique_CDR3_sequences_AA),1,function(x){strsplit(x, split=", ")[[1]]}))))
              for(c in 1:length(cdr3temp)){
                cdr3tempnr<-vector()
                for(g in 1:length(unique(ind.adist[grepCDR3]))){
                  inter<-intersect(grep(unique(ind.adist[grepCDR3])[g],clones.tab.new$samples),grep(gsub("[*]","_",cdr3temp[c]),gsub("[*]","_",clones.tab.new$all_CDR3_sequences_AA)))
                  if(length(inter)==0){
                    cdr3tempnr<-c(cdr3tempnr, paste(unique(ind.adist[grepCDR3])[g], 0, sep=": "))
                  }else{
                    cdr3tempnr<-c(cdr3tempnr, paste(unique(ind.adist[grepCDR3])[g], length(gregexpr(gsub("[*]","_",cdr3temp[c]),gsub("[*]","_",c(clones.tab.new$all_CDR3_sequences_AA[inter])))[[1]])/length(inter), sep=": "))
                  }
                }
                nrseq.temp<-c(nrseq.temp, paste(cdr3temp[g],"(",do.call(paste, c(as.list(as.character(cdr3tempnr)), sep="; ")),")",sep=""))
              }
              nrseq<-do.call(paste, c(as.list(as.character(nrseq.temp)), sep="; "))
              if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                           do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                           uniCDR3length[i],                                                               
                                                           do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist==T))]))), sep="; ")),
                                                           length(unique(CDR3.adist[c(a,which(a.dist==T))])),
                                                           nrseq,
                                                           do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                           clones.uniV.gene[v],
                                                           if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                           
                                                           if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                           if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                           if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                           if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                           if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                           if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                           if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                           if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                           if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                           if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
              } ## end: new entry
              rm(clones.tab.new)
            } ## end if shared CDR3 exists
            rm(a.dist)
          } ## end for adist
          rm(CDR3VJ.index)
          rm(CDR3.adist.temp)
          rm(CDR3.adist)
          rm(ind.adist)
          gc()
        } ## end: use J
      } ## end V genes
    } ## end CDR3 length
      return(sharedclone.temp)
  } ## end dopar

  sharedclone<-do.call(rbind.data.frame, clonelist)
  stopCluster(cl)
  
  if(length(sharedclone)>0){
    if(length(which(as.numeric(sharedclone[,1])==1))>0){
      sharedclone<-sharedclone[which(as.numeric(sharedclone[,1])>1),]
    }
    if(nrow(sharedclone)>0){
      colnames(sharedclone)<-c("number_samples",
                               "samples",
                               "CDR3_length_AA",
                               "shared_CDR3",
                               "number_shared_CDR3",
                               "sequence_count_per_CDR3",
                               "unique_CDR3_sequences_AA",
                               "V_gene",
                               "J_gene", 
                               if(dispD==T){"D_gene"},
                               if(dispCDR3aa==T){"all_CDR3_sequences_AA"},
                               if(dispCDR3nt==T){"all_CDR3_sequences_nt"},
                               if(dispFunctionality.list==T){"Functionality_all_sequences"},
                               if(dispFunctionality.ratio==T){c("Func_productive_sequences","Func_unproductive_sequences")},
                               if(dispJunctionFr.list==T){"Junction_frame_all_sequences"}, 
                               if(dispJunctionFr.ratio==T){c("JF_in_frame","JF_out_of_frame")},
                               if(dispTotalSeq==T){"Total_sequences_nt"})    
      
      sharedclone<-sharedclone[intersect(which(duplicated(sharedclone$samples)==F),which(duplicated(sharedclone$"shared_CDR3")==F)),]
      
    }else{
      sharedclone<-vector()
    }
  }
  return(sharedclone)
}

