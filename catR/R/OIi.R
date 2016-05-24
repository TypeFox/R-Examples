OIi<-function(th,it,x, model=NULL,D=1){
pr<-Pi(th,it,model=model,D=D)
P<-pr$Pi
dP<-pr$dPi
d2P<-pr$d2Pi
if (is.null(model)){
Q<-1-P
res<-(P*Q*dP^2-(x-P)*(P*Q*d2P+dP^2*(P-Q)))/(P^2*Q^2)
}
else{
res<-NULL
for (i in 1:length(x)) res[i]<-dP[i,x[i]+1]^2/P[i,x[i]+1]^2-d2P[i,x[i]+1]/P[i,x[i]+1]
}
return(res)}


