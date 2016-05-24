convmort<-function(value=NULL,fromto=1,type=2,M=NULL){
 if(is.null(value)) stop("value is missing.")
 if(!any(fromto %in% c(1,2))) stop("'fromto' value incorrect.")
 if(!any(type %in% c(1,2))) stop("'type' value incorrect.")
 if(type==2 & is.null(M)) stop("Need M for a Type 2 fishery.")

  if(fromto==1){
      if(type==1){
         mu<-1-exp(-value)
        return(mu)
       }
      if(type==2){
         mu<-(value*(1-exp(-value-M)))/(value+M)
         return(mu)
       }
    }
  if(fromto==2){
      if(type==1){
         F<--log(1-value)
         return(F)
       }
      if(type==2){
            f<-function(x,value,M){
              d<-(x*(1-exp(-x-M)))/(x+M)
              return(abs(value-d))
             } 
         F<-optimize(f,c(0,100),tol=0.0000001,value=value,M=M)[[1]]
         return(F)
       }
   }
}

         

 
