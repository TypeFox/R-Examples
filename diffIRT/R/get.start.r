get.start=function(rt,x,constrain,start,model){
  N=nrow(rt)
  nit=ncol(rt)
  ou=matrix(,N,3)  
  
 pc= apply(x,2,mean,na.rm=T)
 if(T%in%(pc<.5) & model==2) 
  stop("You are trying to fit a Q-diffusion model to response data that 
  includes performance below chance level, P(X=1)<.5. Please provide your 
  own starting values or parameter constraints.\n")
 
  
if(identical(apply((!is.na(x))*1,1,sum),rep(1,nrow(x)))){

        if(constrain[3*nit+1]==0) sap=start[3*nit+1]
        if(constrain[3*nit+2]==0) svp=start[3*nit+2] else svp=.5

        ai.s=vi.s=ter.s=rep(-999,nit)
        
        for(i in 1:nit){
          xx=calcEZ(mean(x[,i],na.rm=T), var(rt[,i],na.rm=T), mean(rt[,i],na.rm=T),N=nrow(x))
          ai.s[i]=1/xx[[2]]
          if(model==1) vi.s[i]=-xx[[1]]
          else vi.s[i]=exp(svp^2/2)/xx[[1]]
          if(xx[[3]]>min(rt[,i],na.rm=T)) ter.s[i]=min(rt[,i],na.rm=T)/2
          else ter.s[i]=xx[[3]]
        }
      
        ter.s[ter.s<0]=1
        strt=strt2=c(ai.s,vi.s,ter.s,sap,svp)
      
        for(i in unique(constrain)) strt2[which(constrain==i)]=mean(strt[which(constrain==i)])
        
        if(strt2[(3*nit+1)]<0) strt2[(3*nit+1)]=0   #this will happen in the case that both vi and sap or svp are fixed 
        if(strt2[(3*nit+1)]<0) strt2[(3*nit+2)]=0
        
        if(model==2) strt2=log(strt2)
        if(model==1) strt2[-((nit+1):(2*nit))]=log(strt2[-((nit+1):(2*nit))])
        return(strt2)
  }


 if(nit>5){
        for(i in 1:N){
          xx=calcEZ(mean(x[i,],na.rm=T), var(rt[i,],na.rm=T), mean(rt[i,],na.rm=T), N=ncol(x))
          ou[i,1]=xx[[1]]
          ou[i,2]=xx[[2]]
          ou[i,3]=xx[[3]]
        }
      
        ou[is.nan(ou)]=NA
      
        vvp=var(ou[,1],na.rm=T)
        vap=var(ou[,2],na.rm=T)
        mvp=mean(ou[,1],na.rm=T)
        map=mean(ou[,2],na.rm=T)
      
        sap=sqrt(log(1+vap/(map^2)))
      
        if(model==1) svp=sqrt(vvp)
        else svp=sqrt(log(1+vvp/(mvp^2)))
      
        if(constrain[3*nit+1]==0) sap=start[3*nit+1]
        if(constrain[3*nit+2]==0) svp=start[3*nit+2]

        ai.s=vi.s=ter.s=rep(-999,nit)
      
        for(i in 1:nit){
          xx=calcEZ(mean(x[,i],na.rm=T), var(rt[,i],na.rm=T), mean(rt[,i],na.rm=T),N=nrow(x))
          ai.s[i]=exp(sap^2/2)/xx[[2]]
          if(model==1) vi.s[i]=-xx[[1]]
          else vi.s[i]=exp(svp^2/2)/xx[[1]]
          if(xx[[3]]>min(rt[,i],na.rm=T)) ter.s[i]=min(rt[,i],na.rm=T)/2
          else ter.s[i]=xx[[3]]
        }
      
        ter.s[ter.s<0]=1
        strt=strt2=c(ai.s,vi.s,ter.s,sap,svp)
      
        for(i in unique(constrain)) strt2[which(constrain==i)]=mean(strt[which(constrain==i)])
        
        if(strt2[(3*nit+1)]<0) strt2[(3*nit+1)]=0   #this will happen in the case that both vi and sap or svp are fixed 
        if(strt2[(3*nit+1)]<0) strt2[(3*nit+2)]=0
        
        if(model==2) strt2=log(strt2)
        if(model==1) strt2[-((nit+1):(2*nit))]=log(strt2[-((nit+1):(2*nit))])
        return(strt2)
  }
   if(nit<6){
  
        ai.s=vi.s=ter.s=rep(-999,nit)
      
        for(i in 1:nit){
          xx=calcEZ(mean(x[,i],na.rm=T), var(rt[,i],na.rm=T), mean(rt[,i],na.rm=T),N=nrow(x))
          ai.s[i]=1/xx[[2]]
          if(model==1) vi.s[i]=-xx[[1]]
          else vi.s[i]=1/xx[[1]]
          if(xx[[3]]>min(rt[,i],na.rm=T)) ter.s[i]=min(rt[,i],na.rm=T)/2
          else ter.s[i]=xx[[3]]
        }
      
        ter.s[ter.s<0]=1
        strt=strt2=c(ai.s,vi.s,ter.s,1,1)
      
        for(i in unique(constrain)) strt2[which(constrain==i)]=mean(strt[which(constrain==i)])
      
        if(model==2) strt2=log(strt2)
        if(model==1) strt2[-((nit+1):(2*nit))]=log(strt2[-((nit+1):(2*nit))])
        return(strt2)
}  
  
}  
  
  
  
  
  


  



