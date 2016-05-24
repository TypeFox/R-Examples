calculateProb <-
function(p, r, x){return(exp(lgamma(r+x)-lgamma(x+1)-lgamma(r)+x*log(p)+r*log(1-p)))}
