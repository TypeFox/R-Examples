summary.cquad <-
function(object, ...){
	
# preliminaries
	out = object
# print output	
	cat("\nCall:\n")
    print(out$call)
    cat("\nLog-likelihood:\n")
    cat(out$lk,"\n")
    tstat = out$coefficients/out$se
    pv = 2*(1-pnorm(abs(tstat)))
    Tab = cbind("est."=out$coefficients,"s.e."=out$se,"t-stat"=tstat,"p-value"=pv)
    print(Tab)
    cat("\n")
    
}
