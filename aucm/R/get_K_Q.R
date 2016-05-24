get.X.diff <- function(x1, ...) UseMethod("get.X.diff") 

get.X.diff.default=function(x1, x2, ...){    
    n1=dim(x1)[1]
    n2=dim(x2)[1]    
    x.diff=x1[rep(1:n1,each=n2), ,drop=FALSE] - x2[rep(1:n2,n1),,drop=FALSE]
}

get.X.diff.formula=function(formula, data, ...){
    
    tmp=model.frame(formula, data)
    X1=model.matrix(formula, tmp[tmp[,1]==1,])[,-1,drop=FALSE]
    X2=model.matrix(formula, tmp[tmp[,1]==0,])[,-1,drop=FALSE]
    n1=nrow(X1)
    n2=nrow(X2)
    X1[rep(1:n1,each=n2), ,drop=FALSE] - X2[rep(1:n2,n1),,drop=FALSE]
}


getQ=function(K, n1, n2, call.C=TRUE, do.pred=FALSE){
    
    if(ncol(K)!=n1+n2) stop("K dimension error: "%+%nrow(K)%+%" by "%+%ncol(K))
    if (!do.pred & nrow(K)!=ncol(K)) stop("K dimension error 2")
    
    if (do.pred) {
        #K is n.pred by n
        
        if (call.C) {
        
            n.pred=nrow(K)
            Q=matrix(0.0,n.pred, n1*n2)
            if (!is.double(K)) K <- as.double(K) #this check is free but we must make
            if (!is.double(Q)) Q <- as.double(Q)
            if (!is.integer(n1)) n1 <- as.integer(n1)
            if (!is.integer(n2)) n2 <- as.integer(n2)
            if (!is.integer(n.pred)) n.pred <- as.integer(n.pred)
            Q = .Call("get_Q_pred", K=K, n1=n1, n2=n2)
            
        } else {
            
            n.pred=nrow(K)            
            i=rep(1:n.pred, n1*n2)
            iq=rep(rep(1:n1,each=n2), each=n.pred)
            jq=rep(rep(1:n2,n1),      each=n.pred)    
            Q=matrix(K[(iq-1)*n.pred+i]-K[(jq+n1-1)*n.pred+i], nrow=n.pred, ncol=n1*n2) # K[i,iq]-K[i,jq+n1]
            
        }
        
    } else {
    
        if (call.C) {
            
            Q=matrix(0.0,n1*n2, n1*n2)
            if (!is.double(K)) K <- as.double(K) #this check is free but we must make
            if (!is.double(K)) Q <- as.double(Q)
            if (!is.integer(n1)) n1 <- as.integer(n1)
            if (!is.integer(n2)) n2 <- as.integer(n2)
            aux=.C("get_Q", K=K, n1=n1, n2=n2, Q=Q, DUP = TRUE, NAOK = FALSE,PACKAGE = "aucm")
            
        } else {
            
            ip=rep(rep(1:n1,each=n2), each=n1*n2)
            jp=rep(rep(1:n2,n1),      each=n1*n2)
            iq=rep(rep(1:n1,each=n2), n1*n2)
            jq=rep(rep(1:n2,n1),      n1*n2)    
            
        #    # another way of forming ip, jp, iq, jq, that is particularly useful in C implementation
        #    p=rep(1:(n1*n2),each=n1*n2)
        #    q=rep(1:(n1*n2),n1*n2)
        #    ip=ceiling(p/n2)
        #    jp=p-(ip-1)*n2
        #    iq=ceiling(q/n2)
        #    jq=q-(iq-1)*n2        
        #    aux2=cbind(ip, jp, iq, jq)            
            
            n=n1+n2
            Q=matrix(K[(ip-1)*n+iq]+K[(jp+n1-1)*n+jq+n1]-K[(ip-1)*n+jq+n1]-K[(jp+n1-1)*n+iq], nrow=n1*n2, ncol=n1*n2) # K[ip,iq]+K[jp+n1,jq+n1]-K[ip,jq+n1]-K[jp+n1,iq]
            # another way of computing Q for linear kernel
            #Q = x.diff %*% diag(1/lambda, ncol(x.diff)) %*%t(x.diff)
            
        }
        
    }
    
    Q
    
}


# ramp function or h
ramp.f=function(eta,s,loss=TRUE){
    if (s==Inf) {
        out=ifelse(eta < 0, 1, 0)
    } else {
        out=ifelse(eta > -.5/s & eta < .5/s, .5-s*eta, 0)
        out[eta <= -.5/s] = 1
    }
    if (loss) drop(out) else drop(1-out)
}
# test ramp.f
# eta=seq(-5,5,length=100); plot(eta, ramp.f(eta,s=1))



h1=function(eta) {
    pmax(1/2-eta,0)
}
h2=function(eta) {
    pmax(-1/2-eta,0)
}
h2.deriv=function(eta){
    -ifelse(eta< -1/2,1,0)
}
