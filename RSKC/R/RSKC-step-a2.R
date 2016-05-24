RSKC.step.a2 <-function(d,gr,ncl,Nouttr,n,f,center){
   # Step (a-2) (trimming) 
   g<-f+1
   dtr<-cbind(d,gr,1:n)
   whoisout<-RSKC.step.a2.trim.Euclidean(dtr,g,ncl,n,center)
   if (Nouttr<1){otr<-n+1}
     else{
      o<-order(whoisout[,1],decreasing=TRUE)[1:Nouttr]
 	  otr<-whoisout[o,2]
    }
 return(sort(otr))
}


RSKC.step.a2.trim.Euclidean<-function(d2,g,ncl,n,mu)
{ # d2 is n by (p+2)
  old<-1;whoisout<-matrix(nrow=n,ncol=2);f<-(g-1);ind<-d2[,g]
  for ( k in 1:ncl)
    {
      dff<-NULL;indW<-NULL;df<-d2[ind==k,-g,drop=FALSE] # df nk by p+1
      dff<-df[,-g,drop=FALSE] # dff nk by p
      Wmu<-mu[k,]
      nk<-nrow(dff);indW<-rowSums(!is.na(dff))
      di<-rowSums(scale(dff,center=Wmu,scale=FALSE)^2,na.rm=TRUE)*f/indW # nk by 1
      new<-nk+old-1
      whoisout[old:new,1:2]<-cbind(di,df[,g])
      old<-new+1
    }  
return(whoisout)
}

