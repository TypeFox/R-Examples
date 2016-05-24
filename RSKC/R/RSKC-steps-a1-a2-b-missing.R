RSKC.a1.a2.b.missing <- function(d,L,ncl,nstart,alpha,n,f,g=f+1,Nout,silent=TRUE)
{
  iteration<-10;
  ## each.trim terminate when ( newss - oldss ) / newss  < e^-8 at i^th step (b) and 
  ## return C from the i^th step(a) and weights from (i-1)^th step (b)
  ## must be Nout=floor(alpha*n)
  oldss<-(-Inf)
  WBSS<-rep(NA,iteration)
  W<-rep(1/sqrt(f),f)
  
  for (i in 1:iteration){	
    ## step a  
    reduced.W<-W[W!=0];sumW<-sum(W)
    w.data<-t(t(d[,W!=0,drop=FALSE])*sqrt(reduced.W)) # reduce the dimention so that the alg is more efficient
    reduced.f<-sum(W!=0)

    A1 <- RSKC.trimkmeans.missing(d=w.data,ncl=ncl,w=reduced.W,trim=alpha,runs=nstart,points=Inf,maxit=10000)

    ## if (Nout!=0){out.a1<-which(A1$classification==(ncl+1))}else{out.a1<-n+1} ## modi Jan 1
    if (Nout!=0){
      out.a1<-A1$oW
    }else{
      out.a1<-n+1
    }

    ##cat("\n\n out.a1 ",out.a1,"length",length(out.a1))
    if (i!=1&oldss!=-Inf)
      {
        Ab4<-RSKC.trimkmeans.missing(d=w.data,ncl=ncl,w=reduced.W,trim=alpha,runs=1,
                                     points=cent.missing(w.data[-save.out.a1,,drop=FALSE],C[-save.out.a1],
                                       ncl,reduced.f),maxit=10000)
        ##if (Nout!=0) out.a1.b4<-which(Ab4$classification==(ncl+1)) else{out.a1.b4<-n+1} ## modi Jan 1
        if (Nout!=0) out.a1.b4 <- Ab4$oW else{out.a1.b4<-n+1}
        
        ##bfsum<-sum(Ab4$disttom[-out.a1.b4]) ## modi Jan 1
        b4sum<-Ab4$WSS
        
        ##afsum<-sum(A1$disttom[-out.a1]) ## modi Jan 1
        afsum<-A1$WSS
        
        
        if (!silent) {p.b4<-round(b4sum,2);p.af<-round(afsum,2)
                      cat("step a: old wwss",p.b4,"new wwss",p.af,"\n")}
        
        if (b4sum<afsum){A1<-Ab4;out.a1<-out.a1.b4}
      }
    ## not necessarily ncl by p mostliely ncl by p' ; p'<<p
    ## w.mu[k,j] is missing if all the cases in C_k in j^th feature are missing   
    w.mu<-A1$means;
    
    ## Assign cluster labels to the outliers
    ## new.C<-A1$classification
    ## new.C can be obtained even if w.mu has missing values
    ## if(Nout!=0) new.C<-class.trimk.missing(w.data,w.mu,new.C,ncl,Nout,reduced.W,sumW) ## modi Jan 1
    new.C<-A1$labels
    
    ## now.mu[k,j] is missing if all the cases in C_k in j^th feature are missing
    ## is.na(now.mu[,W!=0])==is.na(w.mu)
    now.mu<-cent.missing(d[-out.a1,,drop=FALSE],new.C[-out.a1],ncl,f)
    
    if (Nout==0){
      out.a2<-n+1 # skip (a-2) jump to (b)
    }else{
                                        #step a-2 (trimming)
      out.a2<-RSKC.step.a2(d,new.C,ncl,Nout,n,f,now.mu)
      n00<-n-length(union(out.a1,out.a2));
    }
    ## If the reduction of out.a2 makes only one cluster nonempty
    ## repeat step a and a-2 again from the begining.
    C0<-new.C[-c(out.a2,out.a1)]
    if(length(unique(C0))==1){cat("alpha is too large for n \n"); next;}
    save.out.a1<-out.a1;save.out.a2<-out.a2;C<-new.C
    
                                        # step b
                                        # performed wo out.a2 nor out.a1
    B<-RSKC.step.b(L,cbind(d[-c(out.a2,out.a1),],C0),ncl,W,f,n00)
    newss <- B$WBSS 

    if(i!=1&oldss!=-Inf){
      diff<-( newss - oldss ) / newss 
      if(!silent) {p.n<-round(newss,2);p.d<-round(diff,2)
                   cat("step b: wbss",p.n,"wbss dif",p.d,"\n")}
      if( diff < 1e-8 ) break;
    }
    oldss <- newss
    W<-B$We;WBSS[i]<-B$WBSS
  }
  
  return(list(centers=now.mu,oE=save.out.a2,oW=save.out.a1,labels=C,weights=W,WBSS=WBSS[1:(i-1)]))
}

cent.missing<-function(d,C,ncl,f){
  mu<-matrix(NA,ncl,f)
  for ( k in 1 :ncl){
    mu[k,]<-colMeans(d[C==k,,drop=FALSE],na.rm=TRUE)
  }
  return(mu)
}

class.trimk.missing<-function(w.data,mu,trimC,ncl,Nout,w,sumW){
  out.cases<-w.data[trimC==(ncl+1),,drop=FALSE]
  out.WdisC<-WDISC(out.cases,mu,ncl,Nout,w,sumW)
  out.belong.index<-max.col(-out.WdisC)
  trimC[trimC==(ncl+1)]<-out.belong.index
  return(trimC)
}	
