pamr.indeterminate <-  function(prob, mingap=0){
n=nrow(prob)
yhat=rep(NA,n)
for(i in 1:n){

  r=rank(-prob[i,])
if(sum(r==1)==1){
 pr1=prob[i,r==1]
 pr2=prob[i,r==2]
if(pr1-pr2 >= mingap){ yhat[i]=(1:ncol(prob))[r==1]}
}}
yhat=as.factor(dimnames(prob)[[2]][yhat])
return(yhat)
}
