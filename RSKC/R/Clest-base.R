Clest.base <-
function(data,B,alpha,nstart,n,L1=12,maxK,sparse,silent){
# In total, RSKC runs B*(maxK-1)*(1+1) in Cbase one RSKC for learning and another for testing sets
# Within Cbase, CER of partition of learning set by classifer and algorithm is 
# calculated for eack K, B times then median is taken as solution of CER
          ClassError<-matrix(NA,B,(maxK-1));colnames(ClassError)<-c("k=2",3:maxK)
          for (b in 1:B){
              	 if (!silent)cat(" Assessing a random partition", b, "out of",B,"\n")
          #(1) randomly partition the original data to two sets nrow(Lb)=2n/3, nrow(Tb)=n/3
                 Ln<-ceiling(2*n/3);Tn<-n-Ln
                 ind<-sample(1:n,Ln,replace=FALSE)
                 L<-data[ind,];T<-data[-ind,];p<-ncol(L);       
          #(2) apply RSKC to learning set L : return cluster centers and weights(classifiers) for each k
                 cer<-vector(length=(maxK-1))
                 for (ncl in 2:maxK){
                           if (!silent)cat(" K=",ncl,"\n")
                           rL<-RSKC(L,ncl,alpha,L1=L1,nstart=nstart,silent=TRUE)
                           W<-rL$weights
                           sumW<-sum(W)       
          #(3) apply the classifiers to testing set T  
                           trans.L<-reduce.dimention.data(L,W)    
                	   trans.mu<-reduce.dimention.mu(trans.L,W,rL$labels,out=rL$oW,K=ncl,N=Ln)  
                	   trans.T<-reduce.dimention.data(T,W)
                           WdisC<-WDISC(trans.T,trans.mu,ncl,Tn,W[W!=0],sumW)
	                   #(dist btw each obs and center of every cluster)
                           T1<-max.col(-WdisC)  # record the index
          #(4) apply RSKC to testing set T; return classification index   
                           rT<-RSKC(T,ncl,alpha,L1=L1,nstart=nstart,silent=TRUE)
                           T2<-rT$labels
          #(5) compute CER for T
                           ClassError[b,(ncl-1)]<-CER(T1,T2,Tn)
                           }           
                }
cer<-apply(ClassError,2,median,na.rm=TRUE)                
return(list(cer=cer,bootstrap=ClassError)) # return maxK-1 by 1 vector
}

reduce.dimention.data<-function(data,W){
	return(t(t(data[,W!=0,drop=FALSE])*sqrt(W[W!=0]))) 
	}
	
reduce.dimention.mu<-function(reduced.data,W,C,out,K=unique(C),N=nrow(data)){

     if (is.character(out)) out<-N+1
     w.mu<-matrix(NA,K,sum(W!=0))
     C2<-C;C2[out]<--1;C2<-C2[1:N]
     for (k in 1 : K) 
             {
             w.mu[k,]<-colMeans(reduced.data[C2==k,,drop=FALSE],na.rm=TRUE)
             }
	return(w.mu)
	}	
