PR.IARM <-
function(x,y,p,Rmat=diag(k*p),rvec=matrix(0,nrow=k*p)){
x=as.matrix(x); k = ncol(x)
if (k == 1) M = PR2(x,y,p,Rmat,rvec);
if (k > 1) M = PR3(x,y,p,Rmat,rvec);
return(M)}
