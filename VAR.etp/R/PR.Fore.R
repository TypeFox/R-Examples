PR.Fore <-
function(x,y,M,h=10){
k <- ncol(M$AR);
if (k==1) Fore=ARM2.Fore(x,y,M,h)
if (k>1) Fore=ARM3.Fore(x,y,M,h)
return(Fore)}
