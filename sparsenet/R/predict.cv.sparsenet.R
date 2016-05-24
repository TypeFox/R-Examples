predict.cv.sparsenet=function(object,newx,which=c("parms.min","parms.1se"),...){
which=match.arg(which)
switch(which,
       parms.min={lambda=object$parms.min[2];which.gamma=object$which.min[2]},
       parms.1se={lambda=object$parms.1se[2];which.gamma=object$which.1se[2]}
       )
predict(object$sparsenet.fit,newx,s=lambda,which.gamma=which.gamma,...)
}
