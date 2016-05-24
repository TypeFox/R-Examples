mcmc.fmove <- function(xy,t,fdabasis,tpred=t,QQ="CAR2",a=1,b=1,r=1,q=1,n.mcmc=100,num.paths.save=10,sigma.fixed=NA){

    ###################################################
    ## Subroutine
    ###################################################

    rnorm.Q <- function(Q,mu=rep(0,nrow(Q)),X=NA,zero.constraint=FALSE,canon=FALSE,add.eps=FALSE){
        
        n=nrow(Q)
        if(add.eps){
            diag(Q)=diag(Q)*(1+.Machine$double.eps*10)
        }
        L=t(chol(Q))
        
        if(canon==FALSE){
            z=rnorm(n)
            v=solve(t(L),z)
            eta.star=mu+v
        }
        if(canon==TRUE){
            w=solve(L,mu)
            q=solve(t(L),w)
            z=rnorm(n)
            v=solve(t(L),z)
            eta.star=q+v
        }   
        
        qsolve=function(L,B){
            ## solve QX=B where Q=LL^T
            V=solve(L,B)
            X=solve(t(L),V)
            X
        }
        
        if(zero.constraint){
            if(is.na(X)[1]){
                X=Matrix(1,nrow=n,ncol=1)
            }
            V=qsolve(L,X)
            W=t(X)%*%V
            W=forceSymmetric(W)
            L.w=t(chol(W))
            U=qsolve(L.w,t(V))
            U=solve(W)%*%t(V)
            c=t(X)%*%eta.star
            eta=eta.star-t(U)%*%c
        }else{
            eta=eta.star
        }
        
        eta
    }
    
    ###################################################
    ## End Subroutine
    ###################################################
    
    x=xy[,1]
    y=xy[,2]
    xbar=mean(x)
    ybar=mean(y)
    
    ## create the matrix KP that links the data observations to the quasi-continuous path
    KP=eval.basis(t,fdabasis)
    KP=Matrix(KP,sparse=T)

    n=nrow(KP)
    m=ncol(KP)
    
    Psi=eval.basis(tpred,fdabasis)
    Psi=Matrix(Psi,sparse=TRUE)
    
    ## make precision matrix of latent Gaussian path
    m=ncol(KP)
    Q=Diagonal(m-1)
    Q=rBind(Q,0)
    Q=cbind(0,Q)
    Q=-Q-t(Q)
    ones=Matrix(1,nrow=nrow(Q),ncol=1)
    diag(Q)=-as.numeric(Q%*%ones)

    if(QQ=="CAR2"){
        D=Q[-c(1,nrow(Q)),]
        Q=t(D)%*%D
    }

    if(length(dim(QQ))==2){
        Q=QQ
    }
    
    
    KtK=t(KP)%*%KP
    
    s2.save=rep(NA,n.mcmc)
    tau2.save=rep(NA,n.mcmc)

    paths.save=list()
    path.save.idx=round(seq(.5*n.mcmc,n.mcmc,length=num.paths.save))
    idx=1

    s2=40
    tau2=1
    
    ## mcmc
    for(iter in 1:n.mcmc){
        if(iter%%100==0) cat(iter," ")
        ## sample betax and betay
        A=KtK/s2+Q/tau2
        bx=1/s2*t(KP)%*%x
        by=1/s2*t(KP)%*%y
        betax=rnorm.Q(A,bx,canon=TRUE)
        betay=rnorm.Q(A,by,canon=TRUE)
        zx=Psi%*%betax
        zy=Psi%*%betay
        
        ## sample s2
        if(is.na(sigma.fixed)){
            s2=1/rgamma(1,r+n,q+.5*(sum((x-KP%*%betax)^2)+sum((y-KP%*%betay)^2)))
        }else{
            s2=sigma.fixed^2
        }
        
        ## sample tau2
        tau2=1/rgamma(1,a+m,b+.5*as.numeric(t(betax)%*%Q%*%betax+t(betay)%*%Q%*%betay))

        ## save samples
        s2.save[iter]=s2
        tau2.save[iter]=tau2
        if(iter==path.save.idx[idx]){
            paths.save[[idx]]=list(xy=cbind(as.numeric(zx),as.numeric(zy)),t=tpred)
            idx=idx+1
        }
    }
    cat("\n")
    list(s2.save=s2.save,tau2.save=tau2.save,pathlist=paths.save)
}
