ARsensitivity.ci=function(ivmodel, Delta=NULL, conflevel=.95){
if(is.null(ivmodel$x)==TRUE){
print("Refit ivmodel with x=TRUE option in the ivreg function.")
stop()
}
regressors=ivmodel$x$regressors
instruments=ivmodel$x$instruments

# Figure out which columns in regressors and instruments are the same
regressors.same=rep(0,ncol(regressors))
instruments.same=rep(0,ncol(instruments))

for(i in 1:ncol(regressors)){
tempmat=regressors[,i]-instruments
tempsum=apply(abs(tempmat),2,sum)
regressors.same[i]=sum(tempsum==0)>0
if(sum(tempsum==0)>0){
instruments.same[which(tempsum==0)]=i
}
}

if(sum(regressors.same==0)!=1){
print("Input exact one exposure for the sensitivity method.")
stop()
}
if(sum(instruments.same==0)!=1){
print("Input exact one IV for the sensitivity method.")
stop()
}

D=as.matrix(regressors[, which(regressors.same==0)])
X=as.matrix(regressors[, -which(regressors.same==0)])
Z=as.matrix(instruments[, which(instruments.same==0)])
Y=matrix(ivmodel$y, ncol=1)
W=instruments
Zhat.given.Xperp=residuals(lm(Z~X-1))

if(!is.null(Delta)){
  if(is.numeric(Delta) & length(Delta)==2){
    ncp=max(Delta^2)*sum(Zhat.given.Xperp^2)
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }
}else{
ncp=0
}

n=nrow(Y)
k=ncol(X);
q=qf(conflevel, df1=1, df2=n-k-1, ncp=ncp)
cval=q/(n-k-1)

Yhat.given.W=fitted(lm(Y~W-1))
Yhat.given.X=fitted(lm(Y~X-1))
Dhat.given.W=fitted(lm(D~W-1))
Dhat.given.X=fitted(lm(D~X-1))

coef.beta0sq=cval*sum(D^2)-(cval+1)*sum(D*Dhat.given.W)+sum(D*Dhat.given.X)
coef.beta0=-2*cval*sum(D*Y)+2*(cval+1)*sum(Y*Dhat.given.W)-2*sum(Y*Dhat.given.X)
coef.constant=cval*sum(Y^2)-(cval+1)*sum(Y*Yhat.given.W)+sum(Y*Yhat.given.X)
D=coef.beta0^2-4*coef.constant*coef.beta0sq

ci=matrix(NA, ncol=2)
colnames(ci)<-c("lower", "upper")

if(coef.beta0sq==0){
if(coef.beta0>0){
info=c("[",-coef.constant/coef.beta0,",Infinity)")
ci[1,]=c(-coef.constant/coef.beta0, Inf)
type=1
}
if(coef.beta0<0){
info=c("(-Infinity,",-coef.constant/coef.beta0,"]");
ci[1,]=c(-Inf, -coef.constant/coef.beta0)
type=1
}
if(coef.beta0==0){
if(coef.constant>=0){
info="Whole Real Line"
ci[1,]=c(-Inf, Inf)
type=1
}
if(coef.constant<0){
info="Empty Set"
type=2
}
}
}
if(coef.beta0sq!=0){
if(D<0){
if(coef.beta0sq>0){
info="Whole Real Line"
ci[1,]=c(-Inf, Inf)
type=1
}
if(coef.beta0sq<0){
info="Empty Set"
type=2
}
}
if(D>0){
# Roots of quadratic equation
# Roots of quadratic equation
root1=(-coef.beta0+sqrt(D))/(2*coef.beta0sq)
root2=(-coef.beta0-sqrt(D))/(2*coef.beta0sq)
upper.root=max(root1,root2)
lower.root=min(root1,root2)
if(coef.beta0sq<0){
info=paste("[",lower.root,",",upper.root,"]")
ci[1, ]=c(lower.root, upper.root)
type=0
}
if(coef.beta0sq>0){
info= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
ci[1, ]=c(-Inf, lower.root)
ci<-rbind(ci, c(upper.root, Inf))
type=1
}
}
if(D==0){
info="Whole Real Line"
ci[1, ]=c(-Inf, Inf)
type=1
}
}
list(confidence.interval=ci, printinfo=info, ci.type=type)
}


