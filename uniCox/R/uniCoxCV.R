uniCoxCV=function(fit, x,y,status,nfolds=5,folds=NULL){
 this.call <- match.call()

n=length(y)


nlam=length(fit$lamlist)

if(is.null(folds)){
folds=balanced.folds(status,nfolds=nfolds)
}

devcv=ncallcv=matrix(NA,nrow=nfolds,ncol=nlam)
for(i in 1:nfolds){
cat(c("FOLD=",i),fill=T)

ii=folds[[i]]

mx=colMeans(x[-ii,])
xs=scale(x[-ii,],center=mx,scale=F)
v=1/coxvar(t(xs), y[-ii], status[-ii])
s0=quantile(sqrt(v),.5)
xs=scale(x[-ii,],center=F,scale=sqrt(v)+s0)
xval=scale(x[ii,],center=F,scale=sqrt(v)+s0)



junk2=uniCox(xs,y[-ii],status[-ii],nlam=nlam)
val=t(junk2$beta)
 for(j in 1:nrow(val)){
cat(j,fill=T)
 eta.new=rowSums(scale(xval,center=F,scale=1/val[j,]))
 a=coxph(Surv(y[ii],status[ii])~eta.new)
 devcv[i,j]=2*diff(a$loglik)
 ncallcv[i,j]=sum(val[j,]!=0)
}
}
devcvm=colMeans(devcv)
ncallcvm=colMeans(ncallcv)
se.devcvm=sqrt(apply(devcv,2,var)/nfolds)
junk=list(devcvm=devcvm, ncallcvm=ncallcvm, se.devcvm=se.devcvm, devcv=devcv,ncallcv=ncallcv,folds=folds, call=this.call)
class(junk)="uniCoxCVFit"
return(junk)
}






