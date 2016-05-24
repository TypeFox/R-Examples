krm.score.test = function (formula, data, K, regression.type=c("logistic","linear"), verbose=FALSE) {

    t.0=Sys.time()
    
    if (length(regression.type)!=1) stop("regression.type has to be specified")
    regression.type <- match.arg(regression.type)
    
    y=model.frame(formula, data)[,1]
    X=model.matrix(formula, data)
    n=nrow(X)
    
    if (regression.type=="logistic") {    
        
        fit=glm(formula, data, family="binomial")
        
        t.1=Sys.time() # % time from last check point
        
        beta.h=coef(fit)
        mu.h=drop(expit(X %*% beta.h))
        
        # plugin estimator
        D.h = c(mu.h*(1-mu.h)) # in 0.0-4 D.h is a diagonal matrix, it changes to a vector in 0.0-5, but in the comments below, it is still considered diagonal matrix
        extra.kurtosis = (mu.h*(1-mu.h)^4 + (1-mu.h)*mu.h^4) - 3 * D.h**2 
            
        XD.5 <- X * D.h^.5
        . <- eigen(crossprod(XD.5))
        V.beta.h <- crossprod(t(.$vectors) * .$values^(-0.5))  # solve(t(X) %*% D.h %*% X)
        V.eta.h  <- crossprod(t(X %*% .$vectors) * .$values^(-0.5)) # X %*% V.beta.h %*% t(X)
        V.mu.h   <- DXD(D.h,V.eta.h,D.h) # D.h %*% V.eta.h %*% (D.h)
        P.h1=diag(D.h) - V.mu.h
         
        t.2=Sys.time() # % time from last check point
        
        # adjusted estimator
        D.h2 <- D.h + diag(V.mu.h) # D.h2 depends on the previous block
        . <- eigen(crossprod(X * D.h2^.5))
        V.beta.h2 <- crossprod(t(.$vectors) * .$values^(-0.5)) # solve(t(X) %*% D.h2 %*% X)
        V.eta.h2  <- crossprod(t(X %*% .$vectors) * .$values^(-0.5)) # X %*% V.beta.h2 %*% t(X)
        V.mu.h2   <- DXD(D.h2,V.eta.h2,D.h2) # D.h2 %*% V.eta.h2 %*% (D.h2)
        P.h2=diag(D.h2) - V.mu.h2
        
        t.3=Sys.time() # % time from last check point
    
        # terms needed for variance computation
        DV  <- DXD(D.h,V.eta.h,rep(1,n))# D.h %*% V.eta.h
        A.h=diag(n) - DV
        W.h=crossprod(A.h,symprod(K,A.h)) # t(A.h) %*% K %*% A.h
        
        dD.dEta.h = c(mu.h * (1-mu.h) * (1-2*mu.h))
        dD.dBeta.h = dD.dEta.h * X
        KDV <- K %*% DV 
        dTr.dD.h = K + crossprod(DV,KDV) - 2*KDV # K + V.eta.h %*% D.h %*% K %*% D.h %*% V.eta.h - K %*% D.h %*% V.eta.h - V.eta.h %*% D.h %*% K 
        dTr.dBeta.h = diag(dTr.dD.h) %*% dD.dBeta.h
        dTr.dBeta.h_V.beta.h = dTr.dBeta.h %*% V.beta.h
        V.muh1.h = dTr.dBeta.h_V.beta.h %*% t(dTr.dBeta.h)
        
        t.4=Sys.time() # % time from last check point
    
        # estimator of mean
        Ph1K <- P.h1 %*% K
        m1=tr(Ph1K)
        m2=tr(P.h2 %*% K)
        
        # estimator of variance
        V.Q.norm.h = 2*tr(crossprod(Ph1K)) # 2*tr( P.h1 %*% K %*% P.h1 %*% K )
        v2 = varQ (W.h, variance=D.h, extra.kurtosis=extra.kurtosis, do.C=T)
        #v4 = v2 + V.muh1.h - 2 * sum( diag(W.h) * dD.dEta.h * drop(dTr.dBeta.h_V.beta.h %*% t(X)) )
        v4 = v2 + V.muh1.h - 2 * sum( (drop(tcrossprod(dTr.dBeta.h_V.beta.h,X)) * diag(W.h)) * dD.dEta.h )
        if (verbose & v4<0) myprint("v4 is negative")  
        # v4 is V.Qc1.h, could be negative under linear regression
        
        # random quadratic form
        z = y-mu.h
        Q = txSy(z,K,z) # drop(t(y-mu.h) %*% K %*% (y-mu.h))
    #    Qc.h1 = Q - m1
    #    Qc.h2 = Q - m2
            
        # variations of p values
        m=m1; v=V.Q.norm.h; s=v/(2*m); k=2*m^2/v;    p.norm.1=pnorm((Q-m)/sqrt(v)); p.chi.1=pchisq(Q/s, df=k)   
        m=m1; v=v2;         s=v/(2*m); k=2*m^2/v;    p.norm.2=pnorm((Q-m)/sqrt(v)); p.chi.2=pchisq(Q/s, df=k)  # m_I v_I in paper
    #    m=m2; v=v2;         s=v/(2*m); k=2*m^2/v;    p.norm.3=pnorm((Q-m)/sqrt(v)); p.chi.3=pchisq(Q/s, df=k)      
    #    m=m1; v=v4;         s=v/(2*m); k=2*m^2/v;    p.norm.4=pnorm((Q-m)/sqrt(v)); p.chi.4=pchisq(Q/s, df=k)      
        m=m2; v=v4;         s=v/(2*m); k=2*m^2/v;    p.norm.5=pnorm((Q-m)/sqrt(v)); p.chi.5=pchisq(Q/s, df=k)  # m_II v_II in paper
        
        res = c("chiI"=p.chi.2, "chiII"=p.chi.5, "normI"=p.norm.2, "normII"=p.norm.5)    
    
        t.5=Sys.time() # % time from last check point
        if(verbose==2) {
            cat("run time of five blocks in percentage:\n")
            print(round(100*as.numeric(c(t.1-t.0, t.2-t.1, t.3-t.2, t.4-t.3, t.5-t.4))/as.numeric(t.5-t.0)))
        }
        if(verbose==2) myprint(m1, v2, m2, v4)
            
            
    } else if (regression.type=="linear") {
    
        fit=glm(formula, data, family="gaussian")        
        beta.h=coef(fit)
        noise.sd = summary(fit)$dispersion ** .5
        mu.h=drop(X %*% beta.h)
        z = (y-mu.h)/noise.sd
        Q = txSy(z,K,z)  # drop(t(y-mu.h) %*% K %*% (y-mu.h))/sigma2
        # there is a 1e-14 difference bt t(y-mu.h) %*% K %*% (y-mu.h)/noise.sd^2 and (t(y-mu.h)/noise.sd) %*% K %*% (y-mu.h)/noise.sd

        P = diag(n) -  X %*% solve(crossprod(X)) %*% t(X)       
        PK <- P %*% K        
        m1=tr(PK)        
        V.Q.norm.h = 2*tr(crossprod(PK)) # 2*tr( P %*% K %*% P %*% K )
        
        # variations of p values
        m=m1; v=V.Q.norm.h; s=v/(2*m); k=2*m^2/v;    p.norm.1=pnorm((Q-m)/sqrt(v)); p.chi.1=pchisq(Q/s, df=k)   # V.Q.norm.h and v.2 are the same under normal model

        res = c("chiI"=p.chi.1, "chiII"=NA, "normI"=p.norm.1, "normII"=NA)    
        
    }    
    
    res
    
}



# this function is kept here so that when this file is sourced, krm.score.test can run without having to source another file
# mu2: a vector of variance
varQ = function(W, variance, extra.kurtosis, do.C=TRUE) {
    
    if (!is.matrix(W)) W=as.matrix(W)
    if (nrow(W)!=ncol(W)) stop("W is not a square matrix.")
    n=nrow(W)
    if (length(variance)!=n) stop("variance is not of length n.")
    if (length(extra.kurtosis)!=n) stop("extra.kurtosis is not of length n.")
    
    mom=0
    
    if (do.C) {
    
        if (!is.double(W)) W <- as.double(W)    
        aux=.C("var_Q", "_W" = W, "_n"=as.integer(n), "_variance"=as.double(variance), "_extra_kurtosis"=as.double(extra.kurtosis), "_mom"=as.double(mom))
        mom=aux$"_mom"
    
    } else {
        
        for (i in 1:n) {
            mom = mom + W[i,i]^2 * extra.kurtosis[i]
        }
    
        for (i in 1:n) {
            for (j in 1:n) {                        
                mom = mom + 2 * W[i,j]^2 * variance[i] * variance[j] # implements Tr(PKPK) in a different way, see adjusted_score_test_v3.pdf
        }}
    
    }
    
    mom 
    
}
#
## testing
#library(krm) # this has to be run so that var_Q is loaded
#n=1e1
#W=diag(n)
## make it a block diagonal
#for (i in 1:n) {
#    if (i %% 2 ==1) W[i,i+1]=1 else W[i,i-1]=1
#}
#
#varQ (W, variance=rep(1/2,n), extra.kurtosis=rep(1/4,n), do.C=TRUE)
#varQ (W, variance=rep(1/2,n), extra.kurtosis=rep(1/4,n), do.C=FALSE)
