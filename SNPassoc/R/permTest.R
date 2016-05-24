permTest<-function(x,method="minimum",K)
 {
  if (!inherits(x, "WGassociation")) 
    stop("x must be an object of class 'WGassociation'")

  if (is.null(attr(x,"pvalPerm")))
    stop ("\n try again 'scanWGassociation' indicating the number of permutations")
  perms<-attr(x,"pvalPerm")

  type.method<-charmatch(method,c("minimum","rtp"))
  if (is.na(type.method))
    stop ("\n method should be 'minimum' or 'rtp'")

  w<-options(warn=-1)

  if (type.method==1)
   {  
    pmin<-apply(perms,2,min,na.rm=TRUE)

    llhd2 <- function(x,p) { 
      ans <- -sum(log(dbeta(x,p[1],p[2])))
     ans
    }

    pIni<-c(1,100)
    param<-nlm(llhd2,x=pmin,p=pIni)$estimate
    psig<-qbeta(0.05,param[1],param[2])
    valid.SNPs<-length(perms[!is.na(perms[,1]),1])
    discarded.SNPs<-length(perms[is.na(perms[,1]),1])
    ans<-list(pmin=pmin,param=param,psig=psig,nSNPs=c(valid.SNPs,discarded.SNPs))
   }


  if (type.method==2)
   {
     rtpK<-function(x,K)
      {
        x<-sort(x)
        ans<-sum(-log(x[1:K]))
        ans
      }

     perms<-attr(x,"pvalPerm")  
     stat<-rtpK(pvalues(x)[,2],K)
     pval<-apply(perms,2,rtpK,K)
     sig<-mean(pval>=stat)
     valid.SNPs<-length(perms[!is.na(perms[,1]),1])
     discarded.SNPs<-length(perms[is.na(perms[,1]),1])
     ans<-list(rtp=stat,sig=sig,nSNPs=c(valid.SNPs,discarded.SNPs),K=K)
   }

  options(w)
  
  attr(ans,"method")<-type.method
  class(ans)<-"permTest"
  ans
 }

