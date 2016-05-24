"beta.loglik"<-function(X, dam_pos=NULL, sire_pos=NULL, par_pos=NULL, beta=NULL, beta_map=NULL, merge=NULL, mergeN=NULL, nUS=c(0,0)){

   ndamcol<-0
   nsirecol<-0


   if(is.null(X[[1]]$DS)){
     if(is.null(X[[1]]$D)==FALSE){
       ndamcol<-ncol(X[[1]]$D)
     }
     if(is.null(X[[1]]$S)==FALSE){
       nsirecol<-ncol(X[[1]]$S)
     }
     nbeta<-ndamcol+nsirecol
   }else{
     nbeta<-ncol(X[[1]]$DS)
   }

   if(length(beta)==0){
     beta<-matrix(0, nbeta,1)
   }else{
     if(length(beta)!=length(unique(beta_map))){
       stop("beta is wrong size in ped.loglik")
     }else{
       beta<-matrix(beta, length(unique(beta_map)),1)
     }
   }

 llik<-0
 for(i in 1:length(X)){

      beta_tmp<-beta[beta_map]

      if(is.null(merge)==FALSE){
        for(m in 1:length(merge)){
          n1<-mergeN[[i]][,m][1]  
          n2<-mergeN[[i]][,m][2] 
          beta_tmp[m]<-beta_tmp[m]+log(n2/n1)
        }
      }
     
   if(is.null(X[[i]]$DS)){
     if(ndamcol>0){
        if(is.na(dam_pos[i])==FALSE){
          Dpred<-exp(X[[i]]$D%*%beta_tmp[beta_map[1:ndamcol]])
	  l_par<-log(Dpred[dam_pos[i]])
          if(nUS[1]>0){
            Dpred[length(Dpred)]<-Dpred[length(Dpred)]*nUS[1]
          }                   
	  ll_sum<-log(sum(Dpred, na.rm=T))	       
          if(is.na(l_par)==FALSE){
            llik <- llik+l_par-ll_sum
          }
        }
     }
     if(nsirecol>0){    
        if(is.na(sire_pos[i])==FALSE){    
          Spred<-exp(X[[i]]$S%*%beta_tmp[beta_map[ndamcol+c(1:nsirecol)]])
	  l_par<-log(Spred[sire_pos[i]])         
          if(nUS[2]>0){
            Spred[length(Spred)]<-Spred[length(Spred)]*nUS[2]
          }                
	  ll_sum<-log(sum(Spred, na.rm=T))	       
          if(is.na(l_par)==FALSE){
            llik<-llik+l_par-ll_sum
          }
        }
     }
   }else{
     if(is.na(par_pos[i])==FALSE){
       DSpred = exp(X[[i]]$DS%*%beta_tmp)
       l_par<-log(DSpred[par_pos[i]])     
       ll_sum<-log(sum(DSpred, na.rm=T))	       
       if(is.na(l_par)==FALSE){
         llik<-llik+l_par-ll_sum
       }
     }
   }
 }
llik
}





