`clogistCalc` <-
function(n,m,x,set,initb=NA,
    h=0.0001,maxitr=15,epsilon=1e-8,conf.level=0.95){

    alpha<-1-conf.level
    ### first find an initial estimate for beta using unconditional 
    ### logistic regression (if it you can)
    nobs<-length(n)
    uset<-unique(set)
    nsets<- length(uset)
    nbeta<-dim(x)[2]
    betanames<-dimnames(x)[[2]]
    if (  !any(is.na(initb)) )  beta<-initb
    else{
        if (nsets>1 && (nsets+nbeta)<nobs){
            xset<-model.matrix(~as.factor(set))
            xall<-cbind(x,xset)
            glmout<-glm.fit(y=cbind(m,n-m),x=xall,family=binomial())
            beta<-glmout$coefficients[1:nbeta]
        } else if (nsets==1 && (nsets+nbeta)<nobs){
            glmout<-glm.fit(y=cbind(m,n-m),x=x,family=binomial())
            beta<-glmout$coefficients[1:nbeta]
        } else beta<-rep(0,nbeta)
    }
    #print(beta)
    ## center x values within set
    center<-function(x) x - mean(x)
    for (k in 1:nsets){
        kset<-set==uset[k]
        x[kset,]<- apply(as.matrix(x[kset,]),2,center)         
    }
    x<-as.matrix(x)
    ### initialize variables
    loglik0<-0
    beta0<-beta
    ### start iteration
    for (i in 1:maxitr){
        loglik<-0
        finfo<-matrix(0,nbeta,nbeta)
        score<-rep(0,nbeta)
        for (k in 1:nsets){
            kset<-set==uset[k]
            out<-clogistInfo(n[kset],
                m[kset],x[kset,],beta,h)
            finfo<-finfo - out$info
            score<-score+out$score
            loglik<-loglik+out$loglik
        }
        beta<- beta + solve(finfo) %*% score
        betadenom<-beta
        betadenom[abs(beta)<.01]<-.01
        error<- max(  abs( (loglik-loglik0)/loglik ),
              abs( (beta-beta0)/betadenom  ) ) 
        if (error<epsilon) break
        if (i==maxitr) warning(paste("did not converge after ",i,"iterations"))
        beta0<-beta
        loglik0<-loglik
    }
    scoremat<-matrix(rep(NA,nbeta*nsets),nsets,nbeta)
    infoarray<-array(0,c(nsets,nbeta,nbeta))
    loglik<-0  
    for (k in 1:nsets){
        kset<-set==uset[k]
        out<-clogistInfo(n[kset],
            m[kset],x[kset,],beta,h)
        infoarray[k,,]<- - out$info
        scoremat[k,]<-out$score
    }
    names(beta)<-betanames
    output<-list(coefficients=beta,u=scoremat,omega=infoarray)
    output$originalCall<-match.call()
    class(output)<-"clogist"
    output
}

