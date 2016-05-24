pp.regression<-function(x,y=NULL,arg=NULL,residu=NULL,teet=NULL,
h=1,kernel="gauss",M=2,method="poid",argd=NULL,vect=FALSE,seed=1)
{
set.seed(seed)
obsit<-ceiling(runif(M)*n)

d<-dim(x)[2]
n<-dim(x)[1]

if (is.null(residu)){
  residu<-matrix(0,n,M)
  teet<-matrix(0,M,d)
  eval<-matrix(0,n,1)
  for (m in 1:M){
    residu[,m]<-y-eval
    ycur<-y-eval

    obs<-obsit[m]
    argd<-x[obs,]

    theta<-single.index(x,ycur,h=h,method=method,argd=argd)
    teet[m,]<-theta
    xcur<-x%*%theta
    if (!vect){
       estimat<-matrix(0,n,1)
       for (nn in 1:n){
            #arg<-x[nn,]
            #arg<-matrix(arg,d,1)
            #acur<-sum(arg*theta)
            acur<-xcur[nn]
            est<-kernesti.regr(acur,xcur,ycur,h=h,kernel=kernel)
            estimat[nn]<-est
       }
    }
    else{
        estimat<-kernesti.regr(xcur,xcur,ycur,h=h,kernel=kernel,vect=vect)       
    }
    eval<-eval+estimat
  }
}
else eval<-NULL

if (is.null(arg)){ 
  value<-NULL
}
else{
  value<-0
  for (m in 1:M){
     ycur<-matrix(residu[,m],n,1)
     tcur<-matrix(teet[m,],d,1)
     xcur<-x%*%tcur
     arg<-matrix(arg,d,1)
     carg<-sum(arg*tcur)
     w<-kernesti.weights(carg,xcur,h=h,kernel=kernel)
     curvalue<-t(w)%*%ycur
     value<-value+curvalue
  }
}

return(list(eval=eval,residu=residu,teet=teet,value=value))
}



