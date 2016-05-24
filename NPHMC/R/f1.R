f1 <-
function(t,survdist,k,lambda0){
 if(survdist=="exp") {k=1; return(lambda0*k*(lambda0*t)^(k-1)*exp(-(lambda0*t)^k))}
 if(survdist=="weib") {return(lambda0*k*(lambda0*t)^(k-1)*exp(-(lambda0*t)^k))}
}

