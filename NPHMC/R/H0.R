H0 <-
function(t,survdist,k,lambda0){
 if(survdist=="exp") {return(lambda0*t)}
 if(survdist=="weib") {return((lambda0*t)^k)}
}

