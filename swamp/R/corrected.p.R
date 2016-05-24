corrected.p <-
function(feature.assoc,correction="fdr",adjust.permute=T,adjust.rank=T,ties.method="first"){
     if(sum(is.na(c(feature.assoc$observed.p,feature.assoc$permuted.p)))>0){stop("NAs in feat$observed.p or feat$permuted.p")}
     # correction by p.adjust methods
     adjusted.m<-p.adjust(feature.assoc$observed.p,method=correction)
     # correction using permutation set: observed p/expected p 
     if(adjust.permute==T){
     adjusted.s<-rep(NA,(length(feature.assoc$observed.p))) 
     names(adjusted.s)<-names(feature.assoc$observed.p)
     ran<-rank(feature.assoc$observed.p,ties.method=ties.method) ## ties only relevant for when method was AUC
     ranp<-rank(feature.assoc$permuted.p,ties.method=ties.method)
     for(i in 1:length(adjusted.s)){
     adjusted.s[i]<-feature.assoc$observed.p[i]/feature.assoc$permuted.p[which(ranp==ran[i])] # (observed pval/expected pval) for a given rank
     }
     adjusted.s[adjusted.s>1]<-1
     }
     # correction using permutation set: number of false discoveries (permuted p values) for each observed p value as theshold
     if(adjust.rank==T){
     adjusted.t<-rep(NA,(length(feature.assoc$observed.p))) 
     names(adjusted.t)<-names(feature.assoc$observed.p)
     ran<-rank(feature.assoc$observed.p,ties.method=ties.method) ## ties only relevant for when method was AUC
     ranp<-rank(feature.assoc$permuted.p,ties.method=ties.method)
     for(i in 1:length(adjusted.s)){
     adjusted.t[i]<- (sum(feature.assoc$observed.p[i]>feature.assoc$permuted.p))/ran[i]
     }
     adjusted.t[adjusted.t>1]<-1
     }
     return(list(padjust=adjusted.m,adjust.permute=if(adjust.permute==T){adjusted.s},adjust.rank=if(adjust.rank==T){adjusted.t}))
     }

