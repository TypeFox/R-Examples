cent3 <-
function(X,n,m,p,mode){

Y=as.matrix(X)
if (mode==1){
    Y=Cc(Y)
}
if (mode==2){
   Y=permnew(Y,n,m,p)
   Y=Cc(Y)
   Y=permnew(Y,m,p,n)
   Y=permnew(Y,p,n,m)
}
if (mode==3){
   Y=permnew(Y,n,m,p)
   Y=permnew(Y,m,p,n)
   Y=Cc(Y)
   Y=permnew(Y,p,n,m)
}
return(Y)
}
