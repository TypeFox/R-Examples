S0 <-
function(t,pi0,survdist,k,lambda0){
 if(survdist=="exp") {k=1;return(pi0+(1-pi0)*exp(-(lambda0*t)^k))}
 if(survdist=="weib") {return(pi0+(1-pi0)*exp(-(lambda0*t)^k))}
}

