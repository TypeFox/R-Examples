icScoreTest<-function(icFIT,group,scores,alternative="two.sided",tol.svd=10^-8){
    x<-group
    A<-icFIT$A
    n<-dim(A)[[1]]
    model.int<-switch(scores,
        logrank1=1,
        logrank2=2,
        wmw=3)
    if(is.numeric(x)) {
        x <- matrix(x, n, 1)
        if(length(unique(x)) == 2) test <- "2-sample"
        else test <- "correlation"
    } else if(is.character(x) | is.factor(x)) {
        ux <- unique(x)	
        nx <- length(ux)
        xout <- matrix(0, n, nx)
        for(i in 1:nx) {
            xout[x == ux[i], i] <- 1
        }
        x <- xout
        test <- paste(nx, "-sample", sep = "")	
    }
    q <- dim(x)[[2]]
    pf<-icFIT$pf
    gg<-apply(A,1,function(x){sum(x*pf) })
    m<-length(pf)-1
    ## S is survival function, length=m
    ## S[i]=S(s_i)
    S<-1-cumsum(pf[-(m+1)])
    lambda<- pf[-(m+1)]/c(1,S[-m])
    Lambda<-cumsum(lambda)
    ## dpdb.nox[i]=dP(s_{i-1})/dBeta without x
    dpdb.nox<-switch(model.int,
        -S*cumsum( (c(1,S[-m]) - S)/c(1,S[-m]) ),
        S*log(S),
        -S*(1-S))
    calcScores<-function(Ai,mult){
        sum( (Ai[-1] - Ai[-(m+1)])*mult)
    }
    dgdb.nox<-apply(A,1,calcScores,mult=dpdb.nox)
    cc<-dgdb.nox/gg
    mmUpperTri<- matrix(1,m,m)
    mmUpperTri[lower.tri(mmUpperTri)]<-0
    dpdgam<-switch(model.int,
        -(matrix(lambda,m,1) %*% matrix(S,1,m))*mmUpperTri,
        diag(S*log(S)),
        -diag(S*(1-S))
        )
    dgdgam<- matrix(NA,n,m)
    for (j in 1:m){
        dgdgam[,j]<-apply(A,1,calcScores,mult=dpdgam[j,])
    }
    d2pdb2.nox<- switch(model.int,
        S*(Lambda^2 - cumsum(lambda*(1-lambda)) ),
        S*log(S)*(1+log(S)),
        -S*(1-S)*(-1+2*S))       
    d2gdb2.nox<- apply(A,1,calcScores,mult=d2pdb2.nox)
    d2pdbdgam.nox<-switch(model.int,
        (matrix(lambda,m,1) %*% matrix(S*Lambda,1,m) - 
            matrix(lambda*(1-lambda),m,1) %*% matrix(S,1,m))*
            mmUpperTri,
        diag(S*log(S)*(1+log(S))),
        diag(-S*(1-S)*(-1+2*S)) )
    d2gdbdgam.nox<- matrix(NA,n,m)
    for (j in 1:m){
        d2gdbdgam.nox[,j]<-apply(A,1,calcScores,mult=d2pdbdgam.nox[j,])
    }
    d2gdgam2<-array(0,c(n,m,m))
    if (model.int==1){
        for (k in 1:m){
            for (j in 1:m){
                Ik<- Ij<-rep(0,m)
                Ik[k<=(1:m)]<-1
                Ij[j<=(1:m)]<-1
                Ikj<-ifelse(k==j,1,0)
                
                tempMult<- Ik*( Ij*S*lambda[j]*lambda[k] -
                    Ikj*S*lambda[k]*(1-lambda[k]) )
                d2gdgam2[,k,j]<-apply(A,1,calcScores,mult=tempMult)
            }
        }
    } else if (model.int==2){
        for (u in 1:m){
                d2gdgam2[,u,u]<-(A[,u+1] - A[,u])*S[u]*log(S[u])*(1+log(S[u]))
        }
    } else if (model.int==3){
        for (u in 1:m){
                d2gdgam2[,u,u]<-(A[,u+1] - A[,u])*(-S[u])*(1-S[u])*(-1+2*S[u])
        }
    } else stop("model.int must be 1,2, or 3")
    ## d2L.xxx are second derivatives of the log likelihood with respect to parameters
    ##    d2L.dB2 - wrt beta
    ##    d2L.dgam2 - wrt gamma
    ##    d2L.dBdgam - wrt to beta and gamma
    d2L.dB2<-matrix(0,q,q)
    d2L.dgam2 <- matrix(0, m, m)
    d2L.dBdgam <- matrix(0, q, m)
    U<-rep(0,q)
    for (i in 1:n){
        d2L.dB2<-d2L.dB2+ (1/gg[i]) * (matrix(x[i,],q,1) %*% matrix(x[i,],1,q))*(d2gdb2.nox[i] -
            (1/gg[i])*( dgdb.nox[i]^2 ))
        d2L.dBdgam<-d2L.dBdgam +  (1/gg[i]) * (matrix(x[i,],q,1) %*% matrix(
            d2gdbdgam.nox[i,] - (1/gg[i])*dgdb.nox[i]*dgdgam[i,],1,m))
        d2L.dgam2<-d2L.dgam2 + (1/gg[i])*(d2gdgam2[i,,] - 
            (1/gg[i])* matrix(dgdgam[i,],m,1) %*% matrix(dgdgam[i,],1,m) )
        U <- U + x[i,  ] * cc[i]
    }
    ## V is Fishers information variance estimate
    V <-  - (d2L.dB2 - d2L.dBdgam %*% solve(d2L.dgam2) %*% t(d2L.dBdgam))

    if (q==1 | q==2){ 
        Z<-U[1]/sqrt(V[1]) 
        p.lte<-pnorm(Z)
        p.gte<-1-pnorm(Z)
        p.twosidedAbs<- 1-pchisq(Z^2,1)
        # Note for normal theory p-values, p.twosided=p.twosidedAbs
        p.values<-c(p.twosided=p.twosidedAbs,p.lte=p.lte,p.gte=p.gte,p.twosidedAbs=p.twosidedAbs)
        if (alternative=="less" | alternative=="greater"){
            statistic<-Z
            names(statistic)<-"Z"
            parameter<-NULL
        } else {
            statistic<-Z^2
            names(statistic)<-"Chi Square"
            parameter<-1
            names(parameter)<-"df"
        }
        p.value<-switch(alternative,
            less=p.lte,
            greater=p.gte,
            two.sided=p.twosidedAbs,
            two.sidedAbs=p.twosidedAbs)
    } else {
        ## calculate g-inverse
        svdv <- svd(V)
        index <- (svdv$d > tol.svd)
        ginvV <- svdv$v[, index] %*% ((1/svdv$d[index]) * t(svdv$u[, index]))
        chisq.value <- matrix(U,1,q) %*% ginvV %*% matrix(U,q,1)
        df<-q-1
        p.twosided <- 1 - pchisq(chisq.value, df)
        p.values<-c(p.twosided=p.twosided,p.twosidedAbs=p.twosided)
        p.value<-p.twosided
        statistic<-chisq.value
        names(statistic)<-"Chi Square"
        parameter<-df
        names(parameter)<-"df" 
    }
    out<-list(p.value=p.value,p.values=p.values,statistic=statistic,
       parameter=parameter,
       V=V,
       d2L.dB2=d2L.dB2,
       d2L.dgam2=d2L.dgam2,
       d2L.dBdgam=d2L.dBdgam)
    out
}

#set.seed(1)
#x<-rexp(20,1)
#z<-c(rep(0,10),rep(1,10))
#fit<-icfit(x~1)
#fit$A
#icst<-icScoreTest(fit,z,scores="logrank1",alternative="two.sided",tol.svd=10^-8)
#icst
