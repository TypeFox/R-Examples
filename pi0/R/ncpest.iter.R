nparncpp.iter=function(p,estimates=c("all","compromise","pi0","f1"),iter=2,weights,eps=1e-6,keep.cdf=NULL,...)
{
	estimates=match.arg(estimates)
    stopifnot(iter>=0)

    if (missing(weights)){
      if(iter==0 || is.infinite(iter)){
	    weights=1
	    pi0=compromise=f1=Inf
		k=0
	    repeat{
        	fit=nparncpp(p,weights=weights,keep.cdf=keep.cdf,...);
        	k=k+1
        	if(k>=1000 || (abs(pi0-fit$pi0)<eps && 
        			abs(compromise-fit$compromise)<eps && 
        			abs(f1-fit$f1)<eps))break
        	weights=pmin(9999,pmax(0,1/sqrt(fit$pdf.p(fit$bincenters))))
            pi0=fit$pi0;compromise=fit$compromise;f1=fit$f1
        }
        fit_wt=fit
      }else{
      	    weights=1
            for(k in 1:iter){
                fit=nparncpp(p,weights=weights,keep.cdf=keep.cdf,...);
                weights=pmin(9999,pmax(0,1/sqrt(fit$pdf.p(fit$bincenters))))
            }
            fit_wt=fit
      }
    }else
        fit_wt=nparncpp(p,weights=weights,keep.cdf=keep.cdf,...)

	rslt=switch(EXPR=estimates,
		compromise=c(compromise=fit_wt$compromise),
		pi0=c(pi0=fit_wt$pi0),
		f1=c(f1=fit_wt$f1),
		all=
		c(pi0=fit_wt$pi0,compromise=fit_wt$compromise,f1=fit_wt$f1)
	)
    attr(rslt,'iter')=if(exists('k',inherits=FALSE))k else 1
    attr(rslt,'fit')=fit_wt
    rslt
}

