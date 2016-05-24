"mismatches"<-function(X.list, G, mm.tol=999){

     noff<-length(X.list$X)
     ndam<-c(unlist(lapply(X.list$X,function(x){length(x$dam.id)})))	
     nsire<-c(unlist(lapply(X.list$X,function(x){length(x$sire.id)})))	

     if(is.genotype(G[[1]])){
       offid<-as.numeric(names(X.list$X))-1	                                
       damid<-c(unlist(lapply(X.list$X,function(x){x$dam.id})))-1
       sireid<-c(unlist(lapply(X.list$X,function(x){x$sire.id})))-1

       nind<-length(G[[1]])
       nloci<-length(G)
       G<-GtoC(G)
       mmD<-rep(0,sum(ndam))
       mmS<-rep(0,sum(nsire))

  output<-.C("mismatches",
	as.integer(nind),		
        as.integer(noff),              
        as.integer(ndam),              
        as.integer(nsire),            
        as.integer(nloci),		
        as.integer(offid),           
        as.integer(damid),            
        as.integer(sireid),			
        as.integer(mmD),            
        as.integer(mmS),                
        as.integer(G)) 

        sD<-1

         for(i in 1:noff){
          X.list$X[[i]]$mmD<-output[[9]][(sD-1)+1:(ndam[i])]
          sD<-sD+ndam[i]
          dremove<-X.list$X[[i]]$dam.id[which(X.list$X[[i]]$mmD>mm.tol)]
          dremove<-which(X.list$X[[i]]$restdam.id%in%dremove==TRUE)
          if(length(dremove)>0){
            X.list$X[[i]]$restdam.id<-X.list$X[[i]]$restdam.id[-dremove]
          }
        }

        sD<-1

        for(i in 1:noff){
          X.list$X[[i]]$mmS<-output[[10]][(sD-1)+1:(nsire[i])]
          sD<-sD+nsire[i]
          sremove<-X.list$X[[i]]$sire.id[which(X.list$X[[i]]$mmS>mm.tol)]
          sremove<-which(X.list$X[[i]]$restsire.id%in%sremove==TRUE)
          if(length(sremove)>0){
            X.list$X[[i]]$restsire.id<-X.list$X[[i]]$restsire.id[-sremove]
          }
        }
      }else{
        for(i in 1:noff){
          X.list$X[[i]]$mmD<-rep(0,ndam[i])
          X.list$X[[i]]$mmS<-rep(0,nsire[i])
        }
      }
X.list
}
       
