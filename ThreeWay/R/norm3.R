norm3 <-
function(X,n,m,p,mode){

Y=as.matrix(X)
if (mode==1){
    Y=t(nrm2(t(Y)))
	Y=((m*p)^.5)*Y
}
if (mode==2){
    Y=permnew(Y,n,m,p)
    Y=t(nrm2(t(Y)))
    Y=permnew(Y,m,p,n)
    Y=permnew(Y,p,n,m)
    Y=((n*p)^.5)*Y
}
if (mode==3){
    Y=permnew(Y,n,m,p)
    Y=permnew(Y,m,p,n)
    Y=t(nrm2(t(Y)))
    Y=permnew(Y,p,n,m)
    Y=((n*m)^.5)*Y
}
return(Y)
}
