eapSem<-function(thEst,it,x,model=NULL,D=1,priorDist="norm",priorPar=c(0,1),lower=-4,upper=4,nqp=33){
if (is.null(model)){
L<-function(th,it,x) prod(Pi(th,it,D=D)$Pi^x*(1-Pi(th,it,D=D)$Pi)^(1-x))
g<-function(s){
res<-NULL
for (i in 1:length(s)) res[i]<-switch(priorDist,
norm=(s[i]-thEst)^2*dnorm(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
unif=(s[i]-thEst)^2*dunif(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
Jeffreys=(s[i]-thEst)^2*sqrt(sum(Ii(s[i],it,D=D)$Ii))*L(s[i],it,x))
return(res)}
h<-function(s){
res<-NULL
for (i in 1:length(s)) res[i]<-switch(priorDist,
norm=dnorm(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
unif=dunif(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
Jeffreys=sqrt(sum(Ii(s[i],it,D=D)$Ii))*L(s[i],it,x))
return(res)}
X<-seq(from=lower,to=upper,length=nqp)
Y1<-g(X)
Y2<-h(X)
}
else{
LL<-function(th,it,x,model){
prob<-Pi(th,it,model=model,D=D)$Pi
res<-1
for (i in 1:length(x)) res<-res*prob[i,x[i]+1]
return(res)}
gg<-function(s,model){
res <- NULL
for (i in 1:length(s)) res[i]<-switch(priorDist,
norm=(s[i]-thEst)^2*dnorm(s[i],priorPar[1],priorPar[2])*LL(s[i],it,x,model=model), 
unif=(s[i]-thEst)^2*dunif(s[i],priorPar[1],priorPar[2])*LL(s[i],it,x,model=model), 
Jeffreys=(s[i]-thEst)^2*sqrt(sum(Ii(s[i],it,model=model,D=D)$Ii))*LL(s[i],it,x,model=model))
return(res)
}
hh<-function(s,model){
res <- NULL
for (i in 1:length(s)) res[i]<-switch(priorDist, 
norm=dnorm(s[i],priorPar[1],priorPar[2])*LL(s[i],it,x,model=model), 
unif=dunif(s[i],priorPar[1],priorPar[2])*LL(s[i],it,x,model=model), 
Jeffreys=sqrt(sum(Ii(s[i],it,model=model,D=D)$Ii))*LL(s[i],it,x,model=model))
return(res)
}
X<-seq(from=lower,to=upper,length=nqp)
Y1<-gg(X,model=model)
Y2<-hh(X,model=model)
}
RES<-sqrt(integrate.catR(X,Y1)/integrate.catR(X,Y2))
return(RES)}
