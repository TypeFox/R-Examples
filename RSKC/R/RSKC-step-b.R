RSKC.step.b<-
  function(L,d,ncl,W,f,N00){ # d is indexed data
    a<-RSKC.step.b.BSS_js(d,(f+1),ncl,W,N00)
    a<-pmax(a,0)
    de<-BinarySearch.delta(tss.wss=a,L1=L)
                                        # FINALLY obatin the weights which maximize the weighted MBCSS
    SS<-pmax(a-de,0) # p by 1
    Sn2<-norm2(SS)
    W<-SS/Sn2
    WBSS<-sum(W*a)
                                        # the position of inf corresponds to the
                                        # of itr where the ok element corrumption occur 
    return(list(We=W,WBSS=WBSS))
  }


RSKC.step.b.BSS_js<-
  function(dout,g,ncl,W,n){ 
    a<-matrix(NA,nrow=ncl,ncol=g-1);sumW<-sum(W)
    d<-dout[,-g];p<-g-1;
    
    Tmu<-colMeans(d,na.rm=TRUE);Tterm<-scale(d,center=Tmu,scale=FALSE)^2 # n by p
                                        # unless all entries of a column is NA, this Tmu vector has no missing value 
    
    index<-dout[,g] # p by 1
    for ( k in 1:ncl){ 
      if (sum(index==k)==0){
        ##cat("no obs in sm clusters");
        next;
      }else{
          df<-d[index==k,,drop=FALSE] #nk by p
                                        # replace from Rcenter ?
          Wmu<-colMeans(df,na.rm=TRUE) # p by 1
          
          nk<-nrow(df)
          Wterm<-scale(df,center=Wmu,scale=FALSE)^2 #nk by p
          Ttermk<-Tterm[index==k,,drop=FALSE] # nk by p
          
          T_W<-Ttermk-Wterm # nk by p
          
          adjust <- as.vector(sumW / ( (!is.na(Wterm)) %*% W ))
                                        # nk by 1 contains the scaler multiple of each obs
          adjustTerm<-T_W*adjust #nk by p
                                        # each row is multiplied by corresponding entry of adjust vec
          a[k,]<-colSums(adjustTerm,na.rm=TRUE) #p by 1
        }
    }
    a<-colSums(a,na.rm=TRUE) 
                                        # in case sm clus contains no obs, na.rm=TRUE
    return(a)
  }



BinarySearch.delta <- function(tss.wss,L1){
                                        # L1=l1bound
                                        # tss.wss= tssperfeatre - wssperfeature
  if(norm2(tss.wss)==0 || sum(abs(tss.wss/norm2(tss.wss)))<=L1) return(0)
  lam1 <- 0
  lam2 <- max(abs(tss.wss))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(tss.wss,(lam1+lam2)/2)
    if(sum(abs(su/norm2(su)))<L1){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}
