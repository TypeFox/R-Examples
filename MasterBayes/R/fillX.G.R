fillX.G<-function(X.list, A, G, E1=0.005, E2=0.005, marker.type="MSW"){

       mtype.numeric<-sum(c("MSC", "AFLP", "MSW", "SNP")%in%marker.type*c(1:4))

       noff<-length(X.list$X)
       ndam<-c(unlist(lapply(X.list$X,function(x){length(x$dam.id)})))	
       nsire<-c(unlist(lapply(X.list$X,function(x){length(x$sire.id)})))	
       offid<-as.numeric(names(X.list$X))-1	                                
       damid<-c(unlist(lapply(X.list$X,function(x){x$dam.id})))-1
       sireid<-c(unlist(lapply(X.list$X,function(x){x$sire.id})))-1

       nall<-sapply(A, length)
       maxall<-max(nall)
       nloci<-length(nall)
       nind<-length(X.list$id)
       A<-unlist(A)
       if(any(A<1e-10)){
         A[which(A<1e-10)]<-1e-8
       }
       if(is.null(E2)){
         E1<-0.005
       }

       if(is.null(E2)){
         E2<-0.005
       }
  
       G<-GtoC(G, (marker.type!="MSC" & marker.type!="MSW"))
       X_design_G<-rep(0, sum(ndam*nsire))

output<-.C("fillXG",
        as.integer(nind),       # number of individuals sampled
        as.integer(noff),       # number of non-base (offspring) individuals
        as.integer(ndam),       # number of candidate dams per offspring
        as.integer(nsire),      # number of candidate sires per offspring
        as.integer(nloci),	# number of loci
        as.integer(nall),       # number of alleles per locus
        as.integer(maxall),     # number of alleles at most polymorhic locus
        as.integer(offid),      # offspring id
        as.integer(damid),      # candidate dam id's for each offspring
        as.integer(sireid),	# candidate sire id's for each offspring	
        as.double(X_design_G),  # Mendelian transition probabilities dam and sire sampled			
        as.double(A),	# starting allele frequencies
        as.double(E1),	        # starting values of E1 and E2
        as.double(E2),	       
        as.integer(G),           # starting true genotypes    
        as.integer(mtype.numeric)
)

startG<-1
startD<-1
startS<-1

for(i in 1:noff){
X.list$X[[i]]$G<-as.matrix(output[[11]][(startG-1)+1:(ndam[i]*nsire[i])])
startG<-startG+ndam[i]*nsire[i]
}
X.list
}
