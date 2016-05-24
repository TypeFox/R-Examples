



ttest.ppls<-function(ppls.object,ncomp=ppls.object$ncomp.opt,index.lambda=ppls.object$index.lambda){
     jack.object=jack.ppls(ppls.object,ncomp=ncomp,index.lambda=index.lambda)

     my.mean<-coef(jack.object)
     my.sd<-sqrt(diag(vcov(jack.object)))
     my.df=jack.object$k-1
     tvalues <- my.mean/my.sd
    pvalues <- 2 * pt(abs(tvalues), df = my.df, lower.tail = FALSE)
             return(list(pvalues=pvalues,tvalues=tvalues))

}