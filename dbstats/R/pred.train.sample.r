

 ###########################
 #### pred.train.sample ####
 ###########################

 ## Description: fit one dblm model for each individual, using the h bandwidth
 ##     - dist1 is used to find the weights
 ##     - dist2 is used to compute the dblm
 ##
 ##     call kernel.number to achieve the weights S
 ##     call dblm to achieve the generalized linear model
 ##
 ##   returned:
 ##     - Shat: Somoothing hat matrix
 ##     - yhat: fitted values
 ## 
 
 
  
  pred.train.sample<-function(y,dist1,dist2,n,h,h.knn,kind.of.kernel,weights,
            rel.gvar,eff.rank){

     if (is.null(h.knn)){
      hi<-rep(h,n)
     }else{
      hi<-pmax(h.knn+1e-10,h) 
     }

     # initialize yhat, weights S and  Smoothing hat matrix Shat
     yhat <- rep(0,n)
     S<- matrix(0,n,n)
     Shat <- matrix(0,n,n)
     
     if (is.null(eff.rank)) method="rel.gvar"
     else method="eff.rank" 

     # find the fitted values of each individual with dblm model 
     # (effective rank such that the model is estimate with the 95% variability 
     # of the data)
      for (i in 1:n){
       S[i,] <- kernel.number(dist1[i,]^.5/hi[i],j=kind.of.kernel)  
                                   
       iid<-which(S[i,]>0)
       weights_aux<-S[i,iid]
       y_aux<-y[iid]
       dist2_aux<-dist2[iid,iid]
                                          
       class(dist2_aux)<-"D2"
       
       if (!is.null(eff.rank))
         eff.rank_aux<-min(length(y_aux)-1,eff.rank)
       else
         eff.rank_aux<- NULL
         
       aux <- dblm.D2(y=y_aux,D2=dist2_aux,weights=weights[iid]*weights_aux,method=method,
                  rel.gvar=rel.gvar,eff.rank=eff.rank_aux)            
       Hhat<-matrix(0,n,n)
       Hhat[iid,iid]=aux$H
      
       auxshat<- matrix(S[i,]/sum(S[i,]),ncol=n,nrow=n,byrow=TRUE) +
             Hhat%*%(diag(rep(1,n))- matrix(S[i,]/sum(S[i,]),ncol=n,nrow=n,
             byrow=TRUE) ) 
       Shat[i,]<-auxshat[i,]  
       yhat[i] <- aux$fitted.values[which(i==iid)]
     }
     
     # return Shat matrix and yhat values
     return(list(Shat=Shat,yhat=yhat))
  }