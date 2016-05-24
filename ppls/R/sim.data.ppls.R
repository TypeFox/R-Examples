'sim.data.ppls'<-
function(ntrain,ntest,stnr,p,a=NULL,b=NULL){
    if (is.null(a)) a=runif(5,-1,1)
    if (is.null(b)) b=runif(5,-1,1)
    p=max(p,5)
    Xtrain=matrix(runif(p*ntrain,-1,1),ncol=p)
    Xtest=matrix(runif(p*ntest,-1,1),ncol=p)
    fxtrain=0
    fxtest=0
    for (j in 1:5) {
        fxtrain=fxtrain + a[j]*Xtrain[,j] + sin(6*b[j]*Xtrain[,j])
        fxtest=fxtest + a[j]*Xtest[,j] + b[j]*sin(6*Xtest[,j])
        }
    sigma=sqrt(var(fxtrain)/stnr)
    ytrain=fxtrain + rnorm(ntrain,0,sigma)
    ytest=fxtest + rnorm(ntest,0,sigma)
    return(list(Xtrain=Xtrain,ytrain=ytrain,Xtest=Xtest,ytest=ytest,sigma=sigma,a=a,b=b))
}
