`Sigmastar.exp` <-
function(d,logb,logtheta,dist) {return(diag(exp(logb),d) + exp( - exp(logtheta) * dist))}

