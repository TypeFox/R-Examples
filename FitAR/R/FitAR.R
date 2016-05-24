`FitAR` <-
function(z,p,lag.max="default", ARModel="ARz", ...){
if (ARModel=="ARz")
    out<-FitARz(z,p,lag.max=lag.max, ...)
else
    out <- FitARp(z,p,lag.max=lag.max, ...) 
out      
}

