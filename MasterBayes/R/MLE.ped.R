MLE.ped<-function(X.list, ped=NULL, USdam=FALSE, nUSdam=NULL, USsire=FALSE, nUSsire=NULL, threshold=0, checkP=FALSE){

         if(is.null(ped)){
           ped<-matrix(NA, length(X.list$id), 3)
           ped[,1]<-1:length(X.list$id)
         }else{
           ped<-match(ped, X.list$id)
           ped<-matrix(ped, length(X.list$id), 3)                    
         }

         posterior.prob<-rep(1, length(X.list$id))
         acceptdam<-rep(TRUE, length(X.list$id))
         acceptsire<-rep(TRUE, length(X.list$id))

         if(FALSE%in%is.na(ped[,2])){
           getd<-FALSE
         }else{
           getd<-TRUE
         }
         if(FALSE%in%is.na(ped[,3])){
           gets<-FALSE
         }else{
           gets<-TRUE
         }

        if(length(USdam)==1){
          if(USdam==TRUE){
            USdam<-rep(1, length(X.list$X))
            betaDcat<-1
          }else{
            USdam<-NULL
            betaDcat<-NULL
          }
        }else{
          betaDcat<-unique(USdam)
        }

        if(length(USsire)==1){
          if(USsire==TRUE){
            USsire<-rep(1, length(X.list$X))
            betaScat<-1
          }else{
            USsire<-NULL
            betaScat<-NULL
          }
        }else{
          betaScat<-unique(USsire)
        }

        nbetaD<-length(betaDcat)
        nbetaS<-length(betaScat)


         ndam<-unlist(lapply(X.list$X, function(x){length(x$dam.id)}))
         nsire<-unlist(lapply(X.list$X, function(x){length(x$sire.id)}))

         nrestdam<-unlist(lapply(X.list$X, function(x){length(x$restdam.id)}))
         nrestsire<-unlist(lapply(X.list$X, function(x){length(x$restsire.id)}))

         if(length(X.list$X[[1]]$G)==0){
           stop("X.list$X is missing genetic likelihoods")
         }


         for(i in 1:length(X.list$X)){

           d_cat<-match(USdam[i], betaDcat)
           s_cat<-match(USsire[i], betaScat)

           pos_in_id<-as.numeric(names(X.list$X)[i])
 
           if(gets==TRUE | getd==TRUE){

             X<-t(matrix(X.list$X[[i]]$G, nsire[i], ndam[i]))

             if(nbetaS>0){
                X[,nrestsire[i]]<-X[,nrestsire[i]]*nUSsire[s_cat]
             }
             if(nbetaD>0){
                X[nrestdam[i],]<-X[nrestdam[i],]*nUSdam[d_cat]
             }
           }

           if(checkP){
             offid<-pos_in_id
             while(length(offid)>0){
               gp<-which(X.list$X[[i]]$dam.id%in%offid)
               if(length(gp)>0){
                 X[gp,]<-0
               }
               gp<-which(X.list$X[[i]]$sire.id%in%offid)
               if(length(gp)>0){
                 X[,gp]<-0
               }
               offid<-c(ped[,1][ped[,2]%in%offid], ped[,1][ped[,3]%in%offid])
               if(length(offid)>0){
                 offid<-unique(offid)
               }
             }           
           }
     
           if(getd==TRUE & gets==TRUE){      # if neither parents have starting parameteristion
             MLpar<-which.max(t(X))[1]
             MLdam<-ceiling(MLpar/nsire[i])
             MLsire<-MLpar-((ceiling(MLpar/nsire[i])-1)*nsire[i])
             posterior.prob[pos_in_id]<-t(X)[MLpar]/sum(X)
           }else{
             if(gets==FALSE){
               MLsire<-match(ped[,3][pos_in_id], X.list$id[X.list$X[[i]]$sire.id])
               if(is.na(MLsire)==TRUE & nbetaS>0){MLsire<-nsire[i]}  
               if(getd==TRUE){   
                 MLdam<-which.max(X[,MLsire])
                 posterior.prob[pos_in_id]<-X[,MLsire][MLdam]/sum(X[,MLsire])
               }  
             }
             if(getd==FALSE){
               MLdam<-match(ped[,2][pos_in_id], X.list$id[X.list$X[[i]]$dam.id])
               if(is.na(MLdam)==TRUE & nbetaD>0){MLdam<-ndam[i]}
               if(gets==TRUE){ 
                 posterior.prob[pos_in_id]<-X[,MLdam][MLsire]/sum(X[,MLdam])       
                 MLsire<-which.max(X[,MLdam])
               }
             }
           }
           if(getd){
             if(posterior.prob[pos_in_id]<threshold){
               acceptdam[pos_in_id]<-FALSE
             }else{
               if(MLdam==length(X.list$X[[i]]$dam.id) & nbetaD>0){
                 ped[,2][pos_in_id]<-NA
               }else{
                 ped[,2][pos_in_id]<-X.list$X[[i]]$dam.id[MLdam]
               } 
             }
           }
           if(gets){
             if(posterior.prob[pos_in_id]<threshold){
               acceptsire[pos_in_id]<-FALSE
             }else{
               if(MLsire==length(X.list$X[[i]]$sire.id) & nbetaS>0){
                 ped[,3][pos_in_id]<-NA
               }else{
                 ped[,3][pos_in_id]<-X.list$X[[i]]$sire.id[MLsire]
               } 
             }
           }
         }
         ped<-as.character(X.list$id[ped])
         ped<-matrix(ped, length(X.list$id), 3)
         ped[,2][which(acceptdam==FALSE)]<-NA
         ped[,3][which(acceptsire==FALSE)]<-NA
         return(list(P=ped, prob=posterior.prob))
}   
