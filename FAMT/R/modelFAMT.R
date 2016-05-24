modelFAMT <-
function(data,x=1,test=x[1],nbf=NULL,maxnbfactors=8,min.err=1e-03) {
if (class(data)[1]!="FAMTdata") stop("Class of data should be FAMTdata")
if (!is.null(nbf)) optimnbfactors = nbf
if (is.null(nbf)) optimnbfactors = nbfactors(data,x,test,pvalues=NULL,maxnbfactors,min.err)$optimalnbfactors
pval = raw.pvalues(data,x,test)
if (optimnbfactors==0) {
   adjdata = data
   adjpval = pval
   fa = NULL
}
if (optimnbfactors>0) {
   fa = emfa(data,nbf=optimnbfactors,x=x,test=test,pvalues=NULL,min.err=min.err)
   rdata = residualsFAMT(data,x,test,pvalues=NULL)$residuals
   stdev = apply(rdata,2,sd)
   adjdata = data
   adjdata$expression = sweep(data$expression,1,FUN="/",STATS=stdev)-fa$B%*%t(fa$Factors)
   adjpval = raw.pvalues(adjdata,x,test)
   fa = emfa(data,nbf=optimnbfactors,x=x,test=test,pvalues=adjpval$pval,min.err=min.err)
   adjdata$expression = sweep(data$expression,1,FUN="/",STATS=stdev)-fa$B%*%t(fa$Factors)
   adjpval = raw.pvalues(adjdata,x,test)
}

idcovar=data$idcovar
res = list(adjpval=adjpval$pval,adjtest=adjpval$test,adjdata=adjdata,FA=fa,pval=pval$pval,x=x,test=test,
   nbf=optimnbfactors, idcovar=idcovar)

class(res) = c("FAMTmodel","list")
return(res)
}
