PR.order <-
function(x,y,pmax=10){
x=as.matrix(x); k = ncol(x)
if (k == 1) AR.order = ARM.order(x,y,pmax);
if (k > 1) AR.order = ARM2.order(x,y,pmax);
return(AR.order)}
