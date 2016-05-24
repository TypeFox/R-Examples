pgen<-function(est=NULL,limit=NULL,estSD=0,limSD=0,corr=0,dist=1,comp=1,nreps=10000){
  if(is.null(est)) stop ("est is missing.")
  if(is.null(limit)) stop ("limit is missing.")
  if(!is.numeric(est)) stop ("est is not numeric.")
  if(!is.numeric(limit)) stop ("limit is not numeric.")
  dtype<-NULL;ltype<-NULL
  if(length(est)==1){
    dtype<-1
    if(estSD==0) stop ("estSD must be greater than zero.")
  }
  if(length(est)>1) dtype<-2
  if(length(limit)==1) ltype<-1
  if(length(limit)>1) ltype<-2
  biascor<-function(a){
     se<-sd(a)/sqrt(length(a))
     CF<-(se^2)/2
     a2<-exp(a-CF)
     return(a2)
  }
    minv<-NULL; maxv<-NULL;sublimit<-NULL;d1<-NULL;d2<-NULL
   datar<-NULL;outs<-NULL;tot<-NULL;n<-NULL;n1<-NULL;n2<-NULL
   nest<-NULL;nlim<-NULL;cov<-NULL
  # dtype == 1 : est uses mean and SD
 if(dtype==1){ # single est
    if(ltype==1){ # single limit 
        if(limSD>0){ # limit with error
           # solve for cov
               cov<-corr*estSD*limSD
               datar<-mvrnorm(n=nreps, c(est,limit), matrix(c(estSD^2,cov,cov,limSD^2),nrow=2,ncol=2))          
               if(dist==2){
                 datar[,1]<-biascor(datar[,1])
                 datar[,2]<-biascor(datar[,2])
               } 
             if(comp==1)  outs<-length(datar[datar[,1]<datar[,2],2])/nreps
             if(comp==2)  outs<-length(datar[datar[,1]<=datar[,2],2])/nreps
             if(comp==3)  outs<-length(datar[datar[,1]>datar[,2],2])/nreps
             if(comp==4)  outs<-length(datar[datar[,1]>=datar[,2],2])/nreps
        }
        if(limSD==0){ #limit without error
          if(dist==1){
            d1<-rnorm(nreps,est,estSD)
            sublimit<-limit
          } 
          if(dist==2){
            d1<-rlnorm(nreps,est,estSD)
            sublimit<-rlnorm(1,limit,0) 
          }
          if(comp==1)  outs<-length(d1[d1<sublimit])/nreps
          if(comp==2)  outs<-length(d1[d1<=sublimit])/nreps
          if(comp==3)  outs<-length(d1[d1>sublimit])/nreps
          if(comp==4)  outs<-length(d1[d1>=sublimit])/nreps         
       }  
    }
  if(ltype==2){ # limit = individual observations
      n<-length(limit)
      if(dist==1){
        d1 <-rnorm(n,est,estSD)
      } 
      if(dist==2){
        d1<-rlnorm(n,est,estSD)
      } 
      tot <- 0
      for (i in 1:n){
      if(comp==1)  tot <- tot + length(d1[d1<limit[i]])
      if(comp==2)  tot <- tot + length(d1[d1<=limit[i]])
      if(comp==3)  tot <- tot + length(d1[d1>limit[i]])
      if(comp==4)  tot <- tot + length(d1[d1>=limit[i]])
      }
      outs <- tot/(n*n)      
    }  
  }
  
  if(dtype==2){ #est = individual observations
    if(ltype==1){ # limit = single observation
       n1<-length(est)
       if(limSD>0){ # limit with error
         n2<-n1
         if(dist==1) d2 <- rnorm(n2,limit,limSD)
         if(dist==2) {
           d2 <- rlnorm(n2,limit,limSD)
         }
       }
       if(limSD==0){ # limit with no error
        n2<-1
        d2<-limit
       }
       tot <- 0
       for (i in 1:n2){
        if(comp==1)  tot <- tot + length(est[est<d2[i]])
        if(comp==2)  tot <- tot + length(est[est<=d2[i]])
        if(comp==3)  tot <- tot + length(est[est>d2[i]])
        if(comp==4)  tot <- tot + length(est[est>=d2[i]])
       }
       outs <- tot/(n1*n2)      
    }
    if(ltype==2){ #limit = individual observations
      nest<-length(est)
      nlim<-length(limit)
      tot <- 0
      for (i in 1:nlim){
        if(comp==1)  tot <- tot + length(est[est<limit[i]])
        if(comp==2)  tot <- tot + length(est[est<=limit[i]])
        if(comp==3)  tot <- tot + length(est[est>limit[i]])
        if(comp==4)  tot <- tot + length(est[est>=limit[i]])
      }
      outs <- tot/(nest*nlim)  
    }
  }
  return(outs)
}
