`sawsCalc` <-
function(beta,u,omega,test=diag(p),beta0=matrix(0,p,1),conf.level=0.95,method=c("d3","d5","d1","d2","d4","dm"),
    bound=.75){
    method<-match.arg(method,c("d3","d5","d1","d2","d4","dm"))
    p<-length(beta)
    K<-dim(u)[[1]]

    if (is.vector(test)){
        if (any(test==0) & length(test)==p){
            warning("test was vector, but treated as matrix because it had p elements and some were 0")
            test<-matrix(test,1,p)
        } else {
            if (max(test)>p | min(test)<1) stop("test was a vector and was of incorrect form")
            test<-diag(p)[test,]
        }
    } else {
        if (is.matrix(test)){
            if (dim(test)[[2]]!=p) stop("incorrect number of columns for test")
            if (!is.numeric(test)) stop("test matrix should be numeric") 
       } else stop("test not a vector or matrix")
    }

    r<-dim(test)[[1]]
    vminv<- apply(omega,c(2,3),sum)
    vm<- solve(  vminv )
    ## H is a K X p vector, where Hi=diag(H[i,])
    if (method=="d4" | method=="d5"){
        HCalc <-function(omegai,Vm=vm,b=bound){
            (1-pmin(b,diag( omegai %*% Vm )) )^(-.5) }
        H<-t(apply(omega,1,HCalc))
    } else H<-matrix(1,K,p)
    if (method=="d2" | method=="d4"){
        PsiHat<-PsiHatCalc(u,H)
        df<-dfCalc(PsiHat,u,omega,H,test,vm)
    } else if (method=="d3" | method=="d5"){
        PsiTilde<-PsiTildeCalc(u,H,test,omega,vm,vminv,p,K)
        df<-rep(NA,r)
        for (j in 1:r){
            df[j]<-dfCalc(PsiTilde[j,,,],u,omega,H,test[j,],vm)
        }
    } else if (method=="d1" | method=="dm") {
        df<-rep(Inf,r)
    }
    if (method=="dm"){ V<-vm 
    } else if (method=="d1" | method=="d2" | method=="d3"){ 
        V<-vm %*% (t(u) %*% u) %*% vm
    } else if (method=="d4" | method=="d5"){
         V<-vm %*% (t(H*u) %*% (H*u)) %*% vm
       
    }
    se<-sqrt(diag(test %*% V %*% t(test)))
    Tval<- (test %*% (matrix(beta,p,1) - beta0) )/se
    df[df<1]<-1
    p.value<-  1 - pf( (Tval)^2,1,df)
    conf.int<-matrix(NA,r,2)
    conf.int[,1]<- test %*% (matrix(beta,p,1) - beta0)  - se * qt(1-(1-conf.level)/2,df)
    conf.int[,2]<- test %*% (matrix(beta,p,1) - beta0)  + se * qt(1-(1-conf.level)/2,df)
    attr(conf.int,"conf.level")<-conf.level
    estimate<- test %*% (matrix(beta,p,1) - beta0)  
    if (dim(test)[[1]]==p && dim(test)[[2]]==p && all(test==diag(p))){
        if (is.null(dimnames(test)[[1]])) dimnames(test)<- list(names(beta),names(beta)) 
    } else warning("test not identity, confidence intervals and pvalues use beta0")

    out<-list(method=method,test=test,beta0=beta0,coefficients=estimate,df=df,V=V,se=se,t.val=Tval,p.value=p.value,conf.int=conf.int)
    class(out)<-"saws"
    out
}

