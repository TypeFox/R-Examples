"prunePed"<-function(pedigree, keep, make.base=FALSE){

   pedigree2<-pedigree
   pedigree2[,1]<-as.character(pedigree[,1])
   pedigree2[,2]<-as.character(pedigree[,2])
   pedigree2[,3]<-as.character(pedigree[,3])
   keep<-as.character(keep)

   ind.keep<-keep
   nind<-length(ind.keep)+1
   while(length(ind.keep)!=nind){
     nind<-length(ind.keep)
     ind.keep<-union(na.omit(c(unlist(pedigree2[,2:3][match(ind.keep,pedigree2[,1]),]))), ind.keep)
   }
   pedigree2<-pedigree2[sort(match(ind.keep, pedigree2[,1])),]

   if(make.base){

     if(any(match(pedigree2[,2], pedigree2[,1])>match(pedigree2[,1], pedigree2[,1]), na.rm=T)){stop("dams appearing before their offspring: try orderPed form MasterBayes")}
     if(any(match(pedigree2[,3], pedigree2[,1])>match(pedigree2[,1], pedigree2[,1]), na.rm=T)){stop("sires appearing before their offspring: try orderPed from MasterBayes")}

     phenotyped<-pedigree2[,1]%in%keep
     delete<-rep(FALSE, dim(pedigree2)[1])

     for(i in 1:dim(pedigree2)[1]){
       nlinks<-phenotyped[i]+sum(pedigree2[,2]%in%pedigree2[,1][i])+sum(pedigree2[,3]%in%pedigree2[,1][i])+sum(is.na(pedigree2[i,][2:3])==FALSE)    
       if(nlinks<2 & phenotyped[i]==FALSE){                  
         pedigree2[,2][which(as.character(pedigree2[,2])==as.character(pedigree2[,1][i]))]<-NA
         pedigree2[,3][which(as.character(pedigree2[,3])==as.character(pedigree2[,1][i]))]<-NA
         delete[i]<-TRUE                                                            
       }
     } 
     if(any(delete)){
       pedigree2<-pedigree2[-which(delete),]
     }
   }
   pedigree2
   
}
