tutl <-
function(p){
# p-quantile of standard log-Weibull distribution
if (p>1 | p<0) {warning("Argument out of range\n"); return(0)}
log(-log(1-p))}

