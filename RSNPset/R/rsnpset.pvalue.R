rsnpset.pvalue<-function(result, pval.transform=FALSE, qfun=function(x){qvalue(x)$qvalue}) {
    W <-abs(result[[1]][,"W"])
    rk<-abs(result[[1]][,"rank"])
    m<-abs(result[[1]][,"m"])
    
    if (pval.transform==TRUE & attr(result, 'ret.rank')==FALSE){
      message("\nNote: The ranks of the replication variance matrices were not returned.  The degrees of freedom for the replication p-values are based on the ranks of the observed variance matrices.\n")
    }

    p<-pchisq(W,rk,lower.tail=FALSE)
    Q<-qfun(p )

    B<-attr(result, 'B')
    if(B==0) {
        rpv<-data.frame(W,"rank"=rk,m,p,Q)
    }
    else {
        k<-length(W)
        Wbk<-matrix(NA,k,B)
        for(i in 1:B) {
            Wbk[,i]<-result[[1+i]][,"W"]
        }
        rkbk<-matrix(NA,k,B)
        if(attr(result,'ret.rank')==TRUE) {
            for(i in 1:B) { 
                rkbk[,i]<-result[[1+i]][,"rank"]
            }
         }
         else {
            for(i in 1:B) { 
                rkbk[,i]<-result[[1]][,"rank"]
            }
         }
        pB<-vector()
        if(pval.transform==FALSE) {
            pB<-rowSums(Wbk>=W)/B
            pB<-pmax(pB,1/B)		        
            QB=qfun(pB)
            
            rpv<-data.frame(W,"rank"=rk,m,p,pB,Q,QB)
        }
        else {
            xiWbkRbk<-pchisq(Wbk,rkbk,lower.tail=FALSE)
            pB<-rowSums(xiWbkRbk<p)/B

            zeta_b<-apply(abs(xiWbkRbk),2,min,na.rm=T)
            z<-matrix(zeta_b,nrow=k,ncol=B,byrow=TRUE)
            PB<-rowSums(z<p)/B

            pB<-pmax(pB,1/B)
            PB<-pmax(PB,1/B)       
            
            QB=qfun(pB)

            rpv<-data.frame(W,"rank"=rk,m,p,pB,PB,Q,QB)
        }
    }
    rownames(rpv)<-rownames(result[[1]])
    class(rpv)<-c("RSNPset.pvalue","data.frame")
    attr(rpv,"K")<-attr(result,"KAna")
    attr(rpv,"B")<-attr(result, 'B')
    attr(rpv,"pval.transform")<- pval.transform
    return(rpv)
}
