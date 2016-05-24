"extractA"<-function(G,marker.type="MSW"){

   if(is.genotype(G[[1]])==FALSE & is.genotypeD(G[[1]])==FALSE){
     if("id"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="id")]
     }
     if("categories"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="categories")]
     }
     G<-genotype.list(G,marker.type=marker.type)
   }

  A<-lapply(G, function(x){summary(x)$allele.freq[,"Proportion"][1:nallele(x)]})
  if(marker.type=="AFLP" | is.genotypeD(G[[1]])){
  A<-lapply(A, function(x){
       if(any(x==1)){
         warning("some loci are monomorphic: setting allele frequencies to 0.9999 and 0.0001")
         x[which(x==1)]<-0.9999
         x[which(x==0)]<-0.0001
       }else{
         x<-x
       }
       x
     })
  }
  A
}
