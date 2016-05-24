'jack.ppls'<-function(ppls.object,ncomp=ppls.object$ncomp.opt,index.lambda=ppls.object$index.lambda){
	mydims<-dim(ppls.object$coefficients.jackknife)
	ncomp.cv<-mydims[2]
	k<-mydims[4]
         if (ncomp>ncomp.cv){
            ncomp=ncomp.cv
            cat(paste("ncomp is too large and set to ",ncomp,".\n"))
          }
        if (index.lambda>length(ppls.object$lambda)){
          index.lambda=length(ppls.object$lambda)
          cat(paste("index of lambda is too large and set to ",index.lambda," .\n"))
        }
        index.lambda=ppls.object$index.lambda
        p<-length(ppls.object$coefficients)
        l=length(ppls.object$lambda)
        coefficients=ppls.object$coefficients.jackknife[,ncomp,index.lambda,]
        mean.ppls<-apply(coefficients,1,mean)
        vcov.ppls<-cov(t(coefficients)) *((k-1)^2)/k
	outlist=list(coefficients=mean.ppls,covariance=vcov.ppls,k=k,ncomp=ncomp,index.lambda=index.lambda)
	class(outlist)='mypls'
        return(outlist)

}