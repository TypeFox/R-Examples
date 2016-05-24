anderson.rubin.ci=function(ivmodel,conflevel=.95){
y=ivmodel$y
n=length(y);
if(is.null(ivmodel$x)==TRUE){
print("Refit ivmodel with x=TRUE option in the ivreg function.")
stop()
}
regressors=ivmodel$x$regressors
instruments=ivmodel$x$instruments
# Use notation in Davidson and MacKinnon, 2011
# Figure out based on ivmodel fit from ivreg, the elements in Davidson and MacKinnon 
W=instruments;
# Figure out which columns in regressors and instruments are the same
regressors.same.vec=rep(0,ncol(regressors))
instruments.columns.same=rep(0,ncol(regressors)-1)
count=1
for(i in 1:ncol(regressors)){
tempmat=regressors[,i]-instruments
tempsum=apply(abs(tempmat),2,sum)
regressors.same.vec[i]=sum(tempsum==0)>0
if(sum(tempsum==0)>0){
instruments.columns.same[count]=i
count=count+1
}
}

y2column=which(regressors.same.vec==0)
y2=as.matrix(regressors[,y2column])
Z=as.matrix(regressors[,-y2column])
W2=as.matrix(instruments[,-instruments.columns.same])

l=ncol(W);
k=ncol(Z);
q=qf(conflevel,l-k,n-l)
cval=q*(l-k)/(n-l)
y1=matrix(y,ncol=1)

y2hat.given.W=fitted(lm(y2~W-1))
y2hat.given.Z=fitted(lm(y2~Z))
y1hat.given.W=fitted(lm(y1~W-1))
y1hat.given.Z=fitted(lm(y1~Z))
coef.beta0sq=cval*sum(y2^2)-(cval+1)*sum(y2*y2hat.given.W)+sum(y2*y2hat.given.Z)
coef.beta0=-2*cval*sum(y1*y2)+2*(cval+1)*sum(y1*y2hat.given.W)-2*sum(y1*y2hat.given.Z)
coef.constant=cval*sum(y1^2)-(cval+1)*sum(y1*y1hat.given.W)+sum(y1*y1hat.given.Z)
D=coef.beta0^2-4*coef.constant*coef.beta0sq
if(coef.beta0sq==0){
if(coef.beta0>0){
ci=c("[",-coef.constant/coef.beta0,",Infinity)")
}
if(coef.beta0<0){
ci=c("(-Infinity,",-coef.constant/coef.beta0,"]");
}
if(coef.beta0==0){
if(coef.constant>=0){
ci="Whole Real Line"
}
if(coef.constant<0){
ci="Empty Set"
}
}
}
if(coef.beta0sq!=0){
if(D<0){
if(coef.beta0sq>0){
ci="Whole Real Line"
}
if(coef.beta0sq<0){
ci="Empty Set"
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
ci=paste("[",lower.root,",",upper.root,"]")
}
if(coef.beta0sq>0){
ci= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
}
}
if(D==0){
ci="Whole Real Line"
}
}
list(confidence.interval=ci)
}


