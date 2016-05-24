IWLS_CorrZIP <-
function(Loadings,Correlations,corrModel,YList=NULL,
            XList=NULL,ZZIndep=NULL,indepModel=NULL,SSIndep=NULL,BetaList=NULL,Vstart=NULL,OFFSETList=NULL,LinkList=c("Log"),DDRIndep=NULL,DRgammaIndep=NULL,
            RespDist=c("Normal","Normal"),RandDistIndep=NULL,DDY=NULL,DYgamma=NULL,
            FactDist=NULL,FF=NULL,SSF=NULL,CorrMat=NULL,ZZCorr=NULL,RandDistCorr=NULL,DDRCorr=NULL,DRCorrgamma=NULL,CustomVarMat=NULL,SSC=NULL,
            EstimateOverDisp=c(FALSE,FALSE),LaplaceFixed=c(FALSE,FALSE),EstimateCorrelations=TRUE, EstimateVariances=TRUE,StandardErrors=TRUE,
            Info=FALSE,DEBUG=FALSE,CONV=CONV,DRFgamma=NULL,APMethod="REML"){

            # Compose all the designs first #
##            if (!require(Matrix)) stop("Package matrix not installed")
           
            # From design matrix into the two #

            # Maximize TP and BN under independence #
            # Estimate correlation #
           
    # Funkcja do dodawania zer wedlug wierszy i column #
    
    FillZeros<-function(a,nr,nc,up=FALSE,left=FALSE){
        if (nr<nrow(a)) stop("number of specified rows is less than the matrix rows")
        if (nc<ncol(a)) stop("number of specified columns is less than the matrix columns")
        if (nr-nrow(a)>0) zeromatrow<-matrix(0,nr-nrow(a),ncol(a))
        else zeromatrow<-NULL
        if (nc-ncol(a)>0) zeromatcol<-matrix(0,nrow(a),nc-ncol(a))
        else zeromatcol<-NULL
        if (nr-nrow(a)>0 & nc-ncol(a)>0) zeromatdia<-matrix(0,nr-nrow(a),nc-ncol(a))
        else zeromatdia<-NULL
    
    
        if (!(up|left)) b<-rbind(cbind(a,zeromatcol),cbind(zeromatrow,zeromatdia))
        if (up & !left) b<-rbind(cbind(zeromatrow,zeromatdia),cbind(a,zeromatcol))
        if (left & !up) b<-rbind(cbind(zeromatcol,a),cbind(zeromatdia,zeromatrow))
        if (up & left) b<-rbind(cbind(zeromatdia,zeromatrow),cbind(zeromatcol,a))
        b
    }
    # Funkcja laczy po przekatnej #
    dbind<-function(a,b){
        out1<-cbind(a,matrix(0,nrow(a),ncol(b)))
        out2<-cbind(matrix(0,nrow(b),ncol(a)),b)
        out<-rbind(out1,out2)
        out
    }
   
    # Creator of U the inverse of the link for V #
    LinkR<-function(x,RandDist){
        if (RandDist=="Normal") out<-x
        if (RandDist=="Gamma") out<-exp(x)
        if (RandDist=="IGamma") out<-(-1/x)
        if (RandDist=="Beta") out<-exp(x)/(1+exp(x))
        out
    }
   
    # Random effects W vector creator - this takes as an argument u vector#
    WRVC<-function(x,RandDist){
        if (RandDist=="Normal") out<-rep(1,length(x))
        if (RandDist=="Gamma") out<-x
        if (RandDist=="IGamma") out<-x^2
        if (RandDist=="Beta") out<-x*(1-x)
        out
    }
    # x- vscale y- uscale - computes deviances for the estimation of the lambda paramters #
    DevRand<-function(x,y,RandDist){
        if (RandDist=="Normal") out<-y^2
        if (RandDist=="Gamma") out<-2*(y-x-1)
        if (RandDist=="IGamma") out<-2*(log(y)-x-1)
        if (RandDist=="Beta") out<--log(4*y*(1-y))
        out
    }
    DWRDU<-function(x,RandDist){
        if (RandDist=="Normal") out<-rep(0,length(x))
        if (RandDist=="Gamma") out<-rep(1,length(x))
        if (RandDist=="IGamma") out<-2*x
        if (RandDist=="Beta") out<-1-2*x
        out
    }
    D2WRDU2<-function(x,RandDist){
        if (RandDist=="Normal") out<-rep(0,length(x))
        if (RandDist=="Gamma") out<-rep(0,length(x))
        if (RandDist=="IGamma") out<-rep(2,length(x))
        if (RandDist=="Beta") out<-rep(2,length(x))
        return(out)
    }
    # link of the main distribution part - choice between canonical inverse and logarithm #
    LinkY<-function(mu,Link){
        if (Link=="Inverse")    eta<--(1/mu)
        if (Link=="Log")        eta<-log(mu)
        if (Link=="Identity")   eta<-mu
        if (Link=="Logit")      eta<-log(mu/(1-mu))
        if (Link=="Probit")     eta<-qnorm(mu)
        if (Link=="CLogLog")    eta<-log(-log(1-mu))
       eta
    }
    # Inverse of the link #
    InvLinkY<-function(eta,Link){
        if (Link=="Inverse")    mu<--(1/eta)
        if (Link=="Log")        mu<-exp(eta)
        if (Link=="Identity")   mu<-eta
        if (Link=="Logit")      mu<-exp(eta)/(1+exp(eta))
        if (Link=="Probit")     mu<-pnorm(eta)
        if (Link=="CLogLog")    mu<-1-exp(-exp(eta))
        mu
    }
    # Generation of the weight matrix W# # This now has two arguments first one says what is the distribution assumed the second one what is the link #
    # Also added parameter B for binomial distribution #
    
    # WARNING !!!! _ currently only canonical links !!!!!!!!!!!!!!!! #
    # These functions Wmatgen and dWdmugen should be ammended with B factor for the binomial distibution !!!!!!!!!!!!#
    Wmatgen<-function(mu,B,Link,Dist){
        if (Dist=="Normal")     Vmat<-rep(1,length(mu))
        if (Dist=="Poisson")    Vmat<-mu
        if (Dist=="Binomial")   Vmat<-(B-mu)*(mu/B)         # In binomial models mu=p*B therefore the transformation is used g(mu/B)=eta #
        if (Dist=="Gamma")      Vmat<-mu^2  
        if (Dist!="Binomial")   B<-1                        # This makes sure offset is not used here if distribution is different than binomial #
                                                            # Include B everywhere and set it to one for different then binomial distribution 3
        #if (Link=="Inverse")    Wvec<-(1/Vmat)
        #if (Link=="Log")        Wvec<-(1/Vmat)*(mu^2)
        #if (Link=="Identity")   Wvec<-(1/Vmat)*rep(1,length(mu))
        #if (Link=="Logit")      Wvec<-(1/Vmat)*
        Wmat<-Vmat
        Wmat
    }
    
    # Generation of bfuncv #
    bfuncvgen<-function(Vvec,Dist){
        if (Dist=="Normal") out<-((Vvec^2)/2)
        if (Dist=="Gamma")  out<-exp(Vvec)
        if (Dist=="Beta")   out<-log(1+exp(Vvec))
        if (Dist=="IGamma") out<--log(-Vvec)
        return(out)
    }
    
    
    # Still the problem with B in link and variance function for binomial seems not to be solved !!!!!!! #
    
    
    dWdmugen<-function(mu,B,Link,Dist){
        mu1<-mu/B
        if (Dist=="Normal") {
            Vmat<-rep(1,length(mu))
            dVmatdmu<-rep(0,length(mu))
        }
        if (Dist=="Poisson") {
            Vmat<-mu
            dVmatdmu<-rep(1,length(mu))
        }
        if (Dist=="Binomial") {
            Vmat<-(B-mu)*(mu/B)
            dVmatdmu<-1-2*(mu/B)
        }
        if (Dist=="Gamma") {
            Vmat<-mu^2
            dVmatdmu<-2*mu
        }
        if (Dist!="Binomial") B<-1
        
        if (Link=="Inverse")    {
            detadmu <- 1/(mu^2)
            d2etadmu2 <- -2/(mu^3)
        }    
        if (Link=="Log")        {
            detadmu<-1/mu
            d2etadmu2<--1/(mu^2)
            
        }
        if (Link=="Identity")   {
            detadmu<-rep(1,length(mu))
            d2etadmu2<-rep(0,length(mu))
            
        }    
        if (Link=="Logit")      {
            detadmu<-1/(mu*(1-mu1))
            d2etadmu2<--(1-2*mu1)/((mu*(1-mu1))^2)
        }
        
        dWdmu<--(1/Vmat^2)*dVmatdmu*((1/detadmu)^2)+2*(1/Vmat)*(1/detadmu)*(-1/detadmu^2)*d2etadmu2       
        dWdmu
    }
        
    d2Wdmu2gen<-function(mu,B,Link,Dist){
        mu1<-mu/B
        if (Dist=="Normal") {
            Vmat<-rep(1,length(mu))
            dVmatdmu<-rep(0,length(mu))
            d2Vmatdmu2<-rep(0,length(mu))
        }
        if (Dist=="Poisson") {
            Vmat<-mu
            dVmatdmu<-rep(1,length(mu))
            d2Vmatdmu2<-rep(0,length(mu))
        }
        if (Dist=="Binomial") {
            Vmat<-(B-mu)*(mu/B)
            dVmatdmu<-1-2*(mu/B)
            d2Vmatdmu2<--2*(1/B)
        }
        if (Dist=="Gamma") {
            Vmat<-mu^2
            dVmatdmu<-2*mu
            d2Vmatdmu2<-2
        }
        if (Dist!="Binomial") B<-1
        
        if (Link=="Inverse")    {
            detadmu <- 1/(mu^2)
            d2etadmu2 <- -2/(mu^3)
            d3etadmu3 <- 6/(mu^4)
        }    
        if (Link=="Log")        {
            detadmu<-1/mu
            d2etadmu2<--1/(mu^2)
            d3etadmu3<-2/(mu^3)
            
        }
        if (Link=="Identity")   {
            detadmu<-rep(1,length(mu))
            d2etadmu2<-rep(0,length(mu))
            d3etadmu3<-rep(0,length(mu))
            
        }    
        if (Link=="Logit")      {
            detadmu<-1/(mu*(1-mu1))
            d2etadmu2<--(1-2*mu1)/((mu*(1-mu1))^2)
            d3etadmu3<-((2/B)*((mu*(1-mu1))^2)+2*(1-2*mu1)*mu*(1-2*mu1)*(1-mu1))/(mu*(1-mu1))^4
        }
        
        # Add d2Vmatdmu2 and d3etadmu3 to all the functions #
        d2Wdmu2<-2*(1/Vmat^3)*(dVmatdmu^2)*((1/detadmu)^2)-(1/Vmat^2)*d2Vmatdmu2*((1/detadmu)^2)+2*(1/Vmat^2)*dVmatdmu*((1/detadmu)^3)*(d2etadmu2)-
                    2*(1/Vmat^2)*(dVmatdmu)*(1/detadmu)*(-1/detadmu^2)*d2etadmu2-2*(1/Vmat)*(1/detadmu^2)*d2etadmu2*(-1/detadmu^2)*d2etadmu2-
                    4*(1/Vmat)*(1/detadmu)*(-1/detadmu^3)*(d2etadmu2^2)+2*(1/Vmat)*(1/detadmu)*(-1/detadmu^2)*d3etadmu3 
        return(d2Wdmu2)     
    }
    
    # Copy of a function for canonical links - direct computation much easier #
    #d2Wdmu2gen<-function(mu,B,Link,Dist){
    #    if (Dist=="Normal")     Vmat<-rep(0,length(mu))
    #    if (Dist=="Poisson")    Vmat<-rep(0,length(mu))
    #    if (Dist=="Binomial")   Vmat<--2/B        # In binomial models mu=p*B therefore the transformation is used g(mu/B)=eta #
    #    if (Dist=="Gamma")      Vmat<-rep(2,length(mu))  
    #    if (Dist!="Binomial")   B<-1                        # This makes sure offset is not used here if distribution is different than binomial #
                                                            # Include B everywhere and set it to one for different then binomial distribution #
        #if (Link=="Inverse")    Wvec<-(1/Vmat)
        #if (Link=="Log")        Wvec<-(1/Vmat)*(mu^2)
        #if (Link=="Identity")   Wvec<-(1/Vmat)*rep(1,length(mu))
        #if (Link=="Logit")      Wvec<-(1/Vmat)*
    #    Wmat<-Vmat
    #    Wmat    
    
    
    #}
    
    
    
    # Generation of the derivative dmudeta #
    dmudetagen<-function(mu,B,Link,Dist){
        if (Link=="Inverse")    dmudeta<-mu^2
        if (Link=="Log")        dmudeta<-mu
        if (Link=="Identity")   dmudeta<-rep(1,length(mu))
        if (Link=="Logit")      dmudeta<-(B-mu)*(mu/B)
        dmudeta
    }
    # Generation of the derivative dAdmu (y-mu)/Phi is outside#
    #dAdmugen<-function(mu,Link){
    #    if (Link=="Inverse")    dAdmu<-rep(0,length(mu))
    #    if (Link=="Log")        dAdmu<-(1/mu^2)
    #    dAdmu
    #}
   
    # These functions are with fact index # - lets keep it for now
    # These functins must be modified for large matricies
    
    # These matricies must be reprogrammed e.g. nrand and nfact must be replaced by nrandcor and nrandind
     SolverShort<-function(ISIGMAMvec,Z){
        nr<-nrow(Z)
        nc<-ncol(Z)
        SigmaE<-ISIGMAMvec[1:nr]
        SigmaR<-ISIGMAMvec[(nr+1):(nr+nc)] # This wont be a diagonal anymore in the correlated random effects models
        tempmat<-t(Z*SigmaE)%*%Z+diag(SigmaR)
        Inverse<-solve(tempmat)
        rm(tempmat)
        PP2<-cbind(Z%*%Inverse%*%t(Z),Z%*%Inverse)
        PP2<-rbind(PP2,cbind(Inverse%*%t(Z),Inverse))
        PP2<-t(t(PP2)*ISIGMAMvec)
        DiagPP2<-diag(PP2)
        rm(PP2)
        list(Inverse=Inverse,DiagPP2=DiagPP2)
    }

    SolverLong<-function(ISIGMAMvec,zTot){
        SigmaE<-ISIGMAMvec[1:ntot]
        SigmaR<-ISIGMAMvec[(ntot+1):(ntot+qcum[nrandcor+nrandind+1])]
        if (!exists("INV1")) {INV1<-SolverShort(ISIGMAMvec,Z)$Inverse}
        INV1<-as.matrix(INV1)
        AA<-as.matrix(t(X*SigmaE)%*%X)
        BB<-as.matrix(t(X*SigmaE)%*%Z)
        CC<-as.matrix(t(Z*SigmaE)%*%X)
        AA1<-as.matrix(solve(AA-BB%*%INV1%*%CC))
        BB1<--AA1%*%BB%*%INV1
        CC1<--INV1%*%CC%*%AA1
        DD1<-INV1+INV1%*%CC%*%AA1%*%BB%*%INV1
        Inverse<-rbind(cbind(AA1,BB1),cbind(CC1,DD1))
        DPMAT<-rep(0,ntot+qcum[nrandcor+nrandind+1])
        # If n is large do the iteration over the row index #
        DPMAT[1:ntot]<-diag(X%*%AA1%*%t(X*SigmaE)+Z%*%(CC1)%*%t(X*SigmaE)+X%*%(BB1)%*%t(Z*SigmaE)+Z%*%(DD1)%*%t(Z*SigmaE))
        # For the random part #
        DPMAT[(ntot+1):length(DPMAT)]<-diag(DD1)*SigmaR
        rm(AA);rm(BB);rm(CC);rm(INV1)
        tempmat<-rbind(cbind(X,as.matrix(Z)),cbind(matrix(0,qcum[nrandcor+nrandind+1],ptot),diag(qcum[nrandcor+nrandind+1])))
        HELP1<-Inverse%*%t(tempmat*ISIGMAMvec)
        NewParms<-HELP1%*%zTot
        rm(Inverse);rm(HELP1);rm(tempmat)
        rm(AA1);rm(BB1);rm(CC1);rm(DD1)
        list(NewParms=NewParms,DiagPMAT=DPMAT)
    }   
   
         
    nearPD<-function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
        maxits = 100)
    {
        if (!(is.numeric(M) && is.matrix(M) && identical(M, t(M))))
            stop("Input matrix M must be square and symmetric.\n")
        inorm <- function(x) max(rowSums(abs(x)))
        n <- ncol(M)
        U <- matrix(0, n, n)
        X <- M
        iter <- 0
        converged <- FALSE
        while (iter < maxits && !converged) {
            Y <- X
            T <- Y - U
            e <- eigen(Y, symmetric = TRUE)
            Q <- e$vectors
            d <- e$values
            D <- if (length(d) > 1)
                diag(d)
            else as.matrix(d)
            p <- (d > eig.tol * d[1])
            QQ <- Q[, p, drop = FALSE]
            X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
            U <- X - T
            X <- (X + t(X))/2
            conv <- inorm(Y - X)/inorm(Y)
            iter <- iter + 1
            converged <- conv <= conv.tol
        }
        X <- (X + t(X))/2
        e <- eigen(X, symmetric = TRUE)
        d <- e$values
        Eps <- posd.tol * abs(d[1])
        if (d[n] < Eps) {
            d[d < Eps] <- Eps
            Q <- e$vectors
            o.diag <- diag(X)
            X <- Q %*% (d * t(Q))
            D <- sqrt(pmax(Eps, o.diag)/diag(X))
            X[] <- D * X * rep(D, each = n)
        }
        (X + t(X))/2
    }               
    
    # First derivative of the SigmaA matrix #
    # invSigmaMat and dSigmadlambda is a matrix #
    dDDdranmat <- function(X,Z,dWdmu,Wvec,dvhatdlambda,invSigmaMat,dSigmadlambda,WR,dWRdu){
        uprow <- dWdmu*Wvec*as.vector(Z%*%dvhatdlambda)
        downrow <- -invSigmaMat%*%dSigmadlambda%*%(invSigmaMat*WR)+(invSigmaMat*dWRdu*WR*dvhatdlambda)
        uprow1 <- t(X*uprow)%*%X
        uprow2 <- t(X*uprow)%*%Z
        dorow1 <- t(Z*uprow)%*%X
        dorow2 <- t(Z*uprow)%*%Z+downrow
    
        out <- rbind(cbind(uprow1,uprow2),cbind(dorow1,dorow2))
        return(out)
    }

# Second derivative of the SigmaA matrix with respect to random effects parameters #
    d2DDdranmat2 <- function(X,Z,d2Wdmu2,dWdmu,Wvec,dvhatdlambda1,dvhatdlambda2,d2vhatdlambda12,invSigmaMat,dSigmadlambda1,dSigmadlambda2,d2Sigmadlambda12,WR,dWRdu,d2WRdu2){
        uprow <- d2Wdmu2*Wvec*as.vector(Z%*%dvhatdlambda1)*as.vector(Z%*%dvhatdlambda2)+dWdmu*dWdmu*as.vector(Z%*%dvhatdlambda1)*as.vector(Z%*%dvhatdlambda2)+
                dWdmu*Wvec*as.vector(Z%*%d2vhatdlambda12)
        downrow <- invSigmaMat%*%dSigmadlambda1%*%invSigmaMat%*%dSigmadlambda2%*%(invSigmaMat*WR)+invSigmaMat%*%dSigmadlambda2%*%invSigmaMat%*%dSigmadlambda1%*%(invSigmaMat*WR)-
                invSigmaMat%*%d2Sigmadlambda12%*%(invSigmaMat*WR)-invSigmaMat%*%dSigmadlambda1%*%(invSigmaMat*dWRdu*WR*dvhatdlambda2)-
                invSigmaMat%*%dSigmadlambda2%*%(invSigmaMat*dWRdu*WR*dvhatdlambda1)+(invSigmaMat*d2WRdu2*WR*dvhatdlambda1*WR*dvhatdlambda2)+
                (invSigmaMat*dWRdu*dWRdu*WR*dvhatdlambda1*dvhatdlambda2)+(invSigmaMat*dWRdu*WR*d2vhatdlambda12)
    
        uprow1 <- t(X*uprow)%*%X
        uprow2 <- t(X*uprow)%*%Z
        dorow1 <- t(Z*uprow)%*%X
        dorow2 <- t(Z*uprow)%*%Z+downrow
    
        out <- rbind(cbind(uprow1,uprow2),cbind(dorow1,dorow2))
        return(out)

    }
        
    dvhatdranmat <- function(invTT2,invSigmaMat,dSigmadlambda,Psi,Uvec){
        out <- -invTT2%*%(invSigmaMat%*%(dSigmadlambda%*%(invSigmaMat%*%(Psi-Uvec))))
        return(out)
    }

# Second derivative of vhat with respect to random effects varaince components parameters #
d2vhatdranmat2 <- function(invTT2,Z,Phi,dWdmu,Wvec,dvhatdlambda1,dvhatdlambda2,invSigmaMat,dWRdu,WR,dSigmadlambda1,dSigmadlambda2,Psi,Uvec,d2Sigmadlambda12){
    out1 <- (t(Z*as.vector(1/Phi)*dWdmu*Wvec*as.vector(Z%*%dvhatdlambda1))%*%as.vector(Z%*%dvhatdlambda2)) + invSigmaMat%*%(dWRdu*WR*dvhatdlambda1*dvhatdlambda2) - 
            invSigmaMat%*%(dSigmadlambda1%*%(invSigmaMat%*%as.vector(WR*dvhatdlambda2))) - invSigmaMat%*%(dSigmadlambda2%*%(invSigmaMat%*%as.vector(WR*dvhatdlambda1))) - 
            invSigmaMat%*%(dSigmadlambda1%*%(invSigmaMat%*%(dSigmadlambda2%*%(invSigmaMat%*%(Psi-Uvec))))) - invSigmaMat%*%(dSigmadlambda2%*%(invSigmaMat%*%(dSigmadlambda1%*%(invSigmaMat%*%(Psi-Uvec))))) +
            invSigmaMat%*%(d2Sigmadlambda12%*%(invSigmaMat%*%(Psi-Uvec)))
    out <- -invTT2%*%out1
    return(out)
}

# First derivative of h with respect to random effects variance components #

# This function need to be ammended for different random effects using the c(psi,lambda) terms #
dhdranmatInd <- function(Z,y,mu,Phi,dvhatdlambda,invSigmaMat,Psi,Uvec,Vvec,bfuncv,dSigmadlambda,randist){
    out1 <- dvhatdlambda%*%(t(Z*as.vector(1/Phi))%*%(y-mu)+invSigmaMat%*%(Psi-Uvec))+(Psi*Vvec-bfuncv)%*%(invSigmaMat^2)%*%rep(1,length(Vvec))
    # depending on distribution a residual needs to added form c(psi,lambda) #
    return(out1)
}

# Temporary first derivative #
dhdranmatCorr <- function(Z,y,mu,Phi,dvhatdlambda,invSigmaMat,Psi,Uvec,Vvec,bfuncv,dSigmadlambda){
    out1 <- dvhatdlambda%*%(t(Z*as.vector(1/Phi))%*%(y-mu)+invSigmaMat%*%(Psi-Uvec))+0.5*Vvec%*%invSigmaMat%*%dSigmadlambda%*%invSigmaMat%*%Vvec-
                0.5*sum(diag(invSigmaMat%*%dSigmadlambda))
    return(out1)
}

# Second derivative of h with respect to random effects variance components #
d2hdranmatCorrCorr <- function(Z,y,mu,Phi,d2vhatdlambda12,dvhatdlambda1,dvhatdlambda2,Wvec,invSigmaMat,dSigmadlambda1,dSigmadlambda2,d2Sigmadlambda12,Psi,Uvec,Vvec,bfuncv,WR){
    out <- d2vhatdlambda12%*%t(Z*as.vector(1/Phi))%*%(y-mu)-dvhatdlambda1%*%t(Z*as.vector(Wvec/Phi))%*%Z%*%dvhatdlambda2 +
            d2vhatdlambda12%*%invSigmaMat%*%(Psi-Uvec)-dvhatdlambda1%*%invSigmaMat%*%dSigmadlambda2%*%invSigmaMat%*%(Psi-Uvec)-dvhatdlambda2%*%invSigmaMat%*%dSigmadlambda1%*%invSigmaMat%*%(Psi-Uvec)-
            dvhatdlambda1%*%(invSigmaMat)%*%(WR*dvhatdlambda2)+
            (sqrt(abs(Psi*Vvec-bfuncv))*sign(Psi*Vvec-bfuncv))%*%(invSigmaMat%*%dSigmadlambda1%*%invSigmaMat%*%dSigmadlambda2%*%invSigmaMat+invSigmaMat%*%dSigmadlambda2%*%invSigmaMat%*%dSigmadlambda1%*%invSigmaMat-
                invSigmaMat%*%d2Sigmadlambda12%*%invSigmaMat)%*%(sqrt(abs(Psi*Vvec-bfuncv))*sign(Psi*Vvec-bfuncv))-
                0.5*sum(diag(invSigmaMat%*%d2Sigmadlambda12))+0.5*sum(diag(invSigmaMat%*%dSigmadlambda1%*%invSigmaMat%*%dSigmadlambda2))
    return(out)
}


            #################################################################################
            # We need to include a profile likelihood function for correlations p_v_beta(h) #
            #################################################################################
            
            #####
            ##### Maybe this function needs to update the random effects estimates #
            #####
            
          AdjProfCorrelations<-function(ZF){
                # This function utilized lexical scooping for the higher level parameters #
                Correls<-list(0)
                ZF<-list(ZF)
                for (i in 1:length(ZF)){
                    Correls[[i]]<-(1-exp(2*ZF[[i]]))/(1+exp(2*ZF[[i]]))
                # Create the design matrix - SigmaTot ## do we need to use Cholesky step here - alternatively we could use the multivariate normal distribution #
                }
                # Unfold CorrMat #
                TempCorrMat<-list(0)
                for (i in 1:length(CorrMat)){
                    TempCorrMat[[i]]<-CorrMat[[i]]
                    for (j in 1:length(Correls[[i]])){
                        TempCorrMat[[i]][CorrMat[[i]]==j]<-Correls[[i]][j]
                    }
                    diag(TempCorrMat[[i]])<-1
                    CorrMatOut[[i]]<-TempCorrMat[[i]]
                }
   
                LambdaCorr<-exp(DDRCorr%*%DRCorrgamma)
                # now create correlation matrix for all the random effects #
                     
                for (i in 1:length(CorrMat)){   
                    LambdaLocal<-rep(0,qcorr[i])
                    for (j in 1:qcorr[i]){
                        LambdaLocal[j]<-LambdaCorr[cumindCorrIndex[cumqcorr[i]+j]+1]
                    }
                    SigmaMat[[i]]<-sqrt(LambdaLocal)*t(CorrMatOut[[i]]*sqrt(LambdaLocal))
                }
                # merging the matrices by diagonals #
                for (i in 1:length(CorrMat)){
                    if (i==1) {
                        SigmaTot<-SigmaMat[[1]]%x%diag(lcorr[1])
                        invSigmaTot<-solve(SigmaMat[[1]])%x%diag(lcorr[1])   
                    }
                    else {
                        SigmaTot<-dbind(SigmaTot,SigmaMat[[i]]%x%diag(lcorr[i]))
                        invSigmaTot<-dbind(invSigmaTot,solve(SigmaMat[[i]])%x%diag(lcorr[i]))
                    }
                }
                #############################################
                ##### Add here estimation of V and Beta #####
                #############################################
                invSigmaTotIn<-invSigmaTot
                if (nrandind>0) invSigmaTotIn<-dbind(invSigmaTotIn,diag(ISIGMAMvec[(qcum[cumqcorr[length(CorrMat)+1]+1]+1):qcum[length(qcum)]]))
                    convIn<-10
                    INVTEMP<-solve(t(ZOriginal)%*%(ZOriginal*as.vector(Wvec/Phi))+invSigmaTotIn)%*%t(ZOriginal*as.vector(Wvec/Phi))
                    VTCorrTotIn<-VTCorrTot
                while (convIn>0.01) {
                    OldVTCorrTotIn<-VTCorrTotIn
                    VTCorrTotIn<-solve(t(ZOriginal)%*%(ZOriginal*as.vector(Wvec/Phi))+invSigmaTotIn)%*%t(ZOriginal*as.vector(Wvec/Phi))%*%(zmain-X%*%Beta)
                    convIn<-sum(abs(OldVTCorrTotIn-VTCorrTotIn))
                }                
                VTCorrTot<-VTCorrTotIn
                # At this stage we have new SigmaTot #
                # All what has to be computed here is the multivariate normal distributions of random effects where the correlations occur #
                # and also the logarithm of the determinant where also the correlations are present , the other factors are independent of rho #
                # at least in this function because we do not interate over the rho again to reestimate beta vs drgamma dygamma... #    
##                require(mvtnorm)
                hlikelihood<-0
                eta<-TTOriginal[1:ntot,]%*%as.matrix(c(Beta,VTCorrTot))
                for (i in 1:nModels){
                    mu[(cModelsDims[i]+1):cModelsDims[i+1]]<-B[(cModelsDims[i]+1):cModelsDims[i+1]]*InvLinkY(eta[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]])
                    if (RespDist[i]=="Normal") hlikelihood<-hlikelihood+sum(dnorm(YList[[i]],mu[(cModelsDims[i]+1):cModelsDims[i+1]],sd=sqrt(Phi[(cModelsDims[i]+1):cModelsDims[i+1]]),log=TRUE))
                    if (RespDist[i]=="Poisson") {
## hlikelihood<-hlikelihood+sum(dpois(YList[[i]],mu[(cModelsDims[i]+1):cModelsDims[i+1]],log=TRUE)/Phi[(cModelsDims[i]+1):cModelsDims[i+1]],log=TRUE)
                       temp<-sum((-mu[(cModelsDims[i]+1):cModelsDims[i+1]]+YList[[i]]*log(mu[(cModelsDims[i]+1):cModelsDims[i+1]])-lgamma(YList[[i]]+1))/Phi[(cModelsDims[i]+1):cModelsDims[i+1]])
                       hlikelihood<-hlikelihood+temp
                    }
                    if (RespDist[i]=="Binomial") hlikelihood<-hlikelihood+sum(dbinom(YList[[i]],B[(cModelsDims[i]+1):cModelsDims[i+1]],(mu[(cModelsDims[i]+1):cModelsDims[i+1]]/B[(cModelsDims[i]+1):cModelsDims[i+1]]),log=TRUE))
                    if (RespDist[i]=="Gamma") hlikelihood<-hlikelihood+sum(dgamma(YList[[i]],shape=(1/Phi[(cModelsDims[i]+1):cModelsDims[i+1]]),scale=(mu[(cModelsDims[i]+1):cModelsDims[i+1]]*Phi[(cModelsDims[i]+1):cModelsDims[i+1]])))
                }
                hlikelihood1<-hlikelihood
                for (i in 1:length(CorrMat)){
                    VTemp<-unlist(VTCorrTot)[(qcum[cumqcorr[i]+1]+1):qcum[cumqcorr[i+1]+1]] # Extract empirical bayes corresponding to the correlated effects of CorrMat[[i]]
                    noraneff<-cumqcorr[i+1]-cumqcorr[i]
                    VTemp<-matrix(VTemp,length(VTemp)/noraneff,noraneff)
                    hlikelihood<-hlikelihood+sum(dmvnorm(VTemp,mean=rep(0,noraneff),sigma=SigmaMat[[i]],log=TRUE))      
                }
                if (nrandind>0) {
                    for (i in 1:nrandind) {
                        if (RandDistIndep[i]=="Normal") hlikelihood<-hlikelihood+sum(dnorm(VTCorrTot[(qcum[cumqcorr[length(CorrMat)+1]]+i):qcum[cumqcorr[length(CorrMat)+1]]+i+1],log=TRUE))
                    }
                }
                # REMARK: There was a problem with invSigmaTot - it was set to independent unit matrix which is not true #
                for (i in 1:length(CorrMat)){
                    if (i==1) {
                        SigmaTot<-SigmaMat[[1]]%x%diag(lcorr[1])
                        invSigmaTot<-solve(SigmaMat[[1]])%x%diag(lcorr[1])   
                    }
                    else {
                        SigmaTot<-dbind(SigmaTot,SigmaMat[[i]]%x%diag(lcorr[i]))
                        invSigmaTot<-dbind(invSigmaTot,solve(SigmaMat[[i]])%x%diag(lcorr[i]))
                    }
                }
                hlikelihood2<-hlikelihood
                MIDMAT<-dbind(diag(as.vector(Wvec/Phi)),invSigmaTot)
                if ((qcum[cumqcorr[length(CorrMat)+1]+1]+1)-qcum[length(qcum)]<0) MIDMAT<-dbind(MIDMAT,diag(ISIGMAMvec[(qcum[cumqcorr[length(CorrMat)+1]+1]+1):qcum[length(qcum)]]))
                DD<-t(TTOriginal)%*%MIDMAT%*%TTOriginal
                TTOriginal3<-rbind(ZOriginal,diag(ncol(ZOriginal)))
                DD1<-t(TTOriginal3)%*%MIDMAT%*%TTOriginal3
                hlikelihood3<-hlikelihood-0.5*determinant((DD1/(2*pi)),logarithm=TRUE)$modulus
                hlikelihood4<-hlikelihood-0.5*determinant((DD/(2*pi)),logarithm=TRUE)$modulus
                AdjProfLike<--hlikelihood4
                MIDMAT1<-dbind(diag(as.vector(Wvec/Phi)),0*invSigmaTot)
                BB<-t(TTOriginal)%*%MIDMAT1%*%TTOriginal
                pd<- sum(diag(solve(DD) %*% BB))
                caic<--2*hlikelihood1+2*pd
                res<-list(AdjProfLike,hlikelihood2,hlikelihood3,caic) 
                return(res)
            }
            ##################################################################################
   
    # Create the design system of all #
   
   
    # Check if there are correlated random effects #
    # if (is.null(CorrMat)) EstimCorrelation<-FALSE - this condition is not really necessary
    # Therefore we dont need a newton-raphson to estimate correlations if corrmat does not exist #
     
#    YBN<-(Y==0)
#    YTP<-Y
#    YTP[Y==0]<-NA

#    nBN<-nrow(YBN)
#    nTP<-nrow(YTP)
#    pBN<-ncol(XBN)
#    pTP<-ncol(XTP)
   
    # This is a total X matrix #
   
    # First model goes Binomial, second Truncated Poisson - or extensions to Poisson etc #
   
    # First also we model the correlated random effects
   
#    XTT<-dbind(XBN,XTP)
   
#    nrandTP<-length(ZZTP)                   # number of truncate poisson components
#    nrandBN<-length(ZZBN)                   # number of binomial
#    nrandCR<-length(ZZCorr)                   # number of correlated components
     
#    nrandTT<-nrandTP+nrandBN
   
   nModels<-length(YList)
   ModelsDims<-sapply(YList,nrow)
   cModelsDims<-cumsum(c(0,ModelsDims))
   RandDist<-c(RandDistCorr,RandDistIndep) # This specified the distribution of all random effects 
   ntot<-sum(ModelsDims)
   qcorrels<-sapply(Correlations,length)
   cumqcorrels<-cumsum(c(0,qcorrels))

   # Design for the matrices of random effects ##
   if (!is.null(DDRCorr) & !is.null(DDRIndep)) DDR<-dbind(DDRCorr,DDRIndep)
   if (is.null(DDRCorr) & !is.null(DDRIndep)) DDR<-DDRIndep
   if (!is.null(DDRCorr) & is.null(DDRIndep)) DDR<-DDRCorr
   if (is.null(DDRCorr) & is.null(DDRIndep)) stop("You did not specify any design matrix for random effects!")
   
   # Create the design matrices of X covariates and Y #
   for (i in 1:length(YList)){
        if (i==1) {
            Y<-YList[[1]]
            X<-XList[[1]]
        }    
        else {
            Y<-rbind(Y,YList[[i]])
            X<-dbind(X,XList[[i]])
        }
   }
   
    # Index of the covariates over all the models #
    p<-sapply(XList,ncol)
    ptot<-sum(p) 
    pcum<-cumsum(c(0,p))

    # Create the matrix of covariances - first multiply correlations by standard deviations to create variance covariance matrix#
    if(!is.null(CorrMat)) {

        
        # Unfold CorrMat #
        TempCorrMat<-list(0)
        CorrMatOut<-list(0)
        
        for (i in 1:length(CorrMat)){
            TempCorrMat[[i]]<-CorrMat[[i]]
            for (j in 1:length(Correlations[[i]])){
                TempCorrMat[[i]][CorrMat[[i]]==j]<-Correlations[[i]][j]
            }
            diag(TempCorrMat[[i]])<-1
            CorrMatOut[[i]]<-TempCorrMat[[i]]
        }
   
        LambdaCorr<-exp(DDRCorr%*%DRCorrgamma)
        # now create correlation matrix for all the random effects #
        qcorr<-rep(0,length(CorrMat))
        lcorr<-rep(0,length(CorrMat)) # in this vector we store

        for (i in 1:length(CorrMat)){
            if (i==1) {
                qcorr[1]<-nrow(CorrMat[[1]])   
            }
            else qcorr[i]<-nrow(CorrMat[[i]])
            cumqcorr<-cumsum(c(0,qcorr))
            lcorr[i]<-ncol(ZZCorr[[cumqcorr[i]+1]])
        }
        for (i in 1:length(CorrMat)){
            tempindex<-as.numeric(names(table(corrModel[(cumqcorr[i]+1):cumqcorr[i+1]])))
            if (sum(LaplaceFixed[tempindex]==FALSE)!=length(LaplaceFixed[tempindex]) & sum(LaplaceFixed[tempindex]==TRUE)!=length(LaplaceFixed[tempindex])) 
                stop("You choose for some correlated effect LAPFIX=TRUE while others LAPFIX=FALSE this is not permitted!")
        }
        # create index of individual random effects #
        indCorrIndex<-rep(0,length(ZZCorr))
        for (i in 1:length(ZZCorr)){
            indCorrIndex[i]<-ncol(ZZCorr[[i]])
        }
        cumindCorrIndex<-cumsum(c(0,indCorrIndex))
        
        SigmaMat<-list(0)
        for (i in 1:length(CorrMat)){   
            LambdaLocal<-rep(0,qcorr[i])
            for (j in 1:qcorr[i]){
                LambdaLocal[j]<-LambdaCorr[cumindCorrIndex[cumqcorr[i]+j]+1]
            }
            SigmaMat[[i]]<-sqrt(LambdaLocal)*t(CorrMatOut[[i]]*sqrt(LambdaLocal))
        }
        # merging the matrices by diagonals #
        for (i in 1:length(CorrMat)){
            if (i==1) {
                SigmaTot<-SigmaMat[[1]]%x%diag(lcorr[1])
                invSigmaTot<-solve(SigmaMat[[1]])%x%diag(lcorr[1])   
            }
            else {
                SigmaTot<-dbind(SigmaTot,SigmaMat[[i]]%x%diag(lcorr[i]))
                invSigmaTot<-dbind(invSigmaTot,solve(SigmaMat[[i]])%x%diag(lcorr[i]))
            }
        }
    
    # Matrix SigmaTot is the resulting matrix #
    # We have to make the random effects indendent via Cholesky decomposition #
    # We need a kroneker product of cholesky matrix times the SigmaTot #
    # You can do the cholesky on the total sigma matrix - the problem is the dimension is greater so maybe we loose computational efficiency #
    # DO cholesky on SigmaMat list #
    # We have to make functions to convert the ZZCorr into vectoral design for correlated random effects and back to the diagnoal design according to subject#
    # This will make the computationally more efficient things #
    
        ZZCorrVec<-list(0)
        for (i in 1:length(ZZCorr)){
            ZZCorrVec[[i]]<-ZZCorr[[i]]%*%rep(1,ncol(ZZCorr[[i]]))
        }
    # Determine how many models we have linked by correlation #
   
    
    # Now we modify the design matrix via cholesky decompositions #
    # All these steps need to be reprogramed using matrix package although in the final product there might be not so many zeros sometime #
        ZZCorrUpd<-list(0)
        DiagDesign<-list(0)
        CholeskyMatrices<-list(0)
        ZZShort<-list(0)
        for (i in 1:length(CorrMat)){
            itchol<-t(chol(SigmaMat[[i]]))  # This is actually cholesky decomposition instead of inverse, there was before inverse which was wrong
            CholeskyMatrices[[i]]<-itchol
            currentindex<-seq(cumqcorr[i]+1,cumqcorr[i+1])
            lengthcorrModelCur<-length(RespDist)
            valuescorrModelCur<-as.numeric(names(table(corrModel)))
            ZZCorrTemp<-rep(list(0),lengthcorrModelCur)
            for (j in 1:lengthcorrModelCur){
                for (k in currentindex) {
                    if (ZZCorrTemp[[j]][1]==0 & length(ZZCorrTemp[[j]])==1){
                        if (corrModel[k]==valuescorrModelCur[j]) ZZCorrTemp[[j]]<-ZZCorrVec[[k]] # If the observation belongs to this model than make it #
                        else ZZCorrTemp[[j]]<-rep(0,ModelsDims[j])
                    }
                    else {
                        if (corrModel[k]==valuescorrModelCur[j]) ZZCorrTemp[[j]]<-cbind(ZZCorrTemp[[j]],ZZCorrVec[[k]])
                        else ZZCorrTemp[[j]]<-cbind(ZZCorrTemp[[j]],rep(0,ModelsDims[j]))           
                    }
                }
            }
        # Binding it all together #
            for (j in 1:length(ZZCorrTemp)){
                if (j==1) {
                    ZZCorrTempTot<-ZZCorrTemp[[j]]
                    nrowtot<-nrow(ZZCorrTemp[[j]])           
                }
                else {
                    ZZCorrTempTot<-rbind(ZZCorrTempTot,ZZCorrTemp[[j]])
                    nrowtot<-c(nrowtot,nrow(ZZCorrTemp[[j]]))
                }
            }
            cnrowtot<-cumsum(c(0,nrowtot))
        # Now we use cholesky transform on the design matrix #
            
            ZZCorrTempTotUpd<-ZZCorrTempTot%*%itchol
            ZZShort[[i]]<-ZZCorrTempTot
        # ZZCorrTempTotUpd is the new design matrix for the joint model from the correlated part  #
        # This design matrix is in the short form (vector form) - we need to expand it to the diagnoal form # 
        # Now we need to take into accout to which model we should link the SSC #
        
        
        # Expansion to diagnoal here we need matrix package already #
            
            for (j in currentindex){
                DiagDesign[[j]]<-model.matrix(~SSC[[j]]-1)
                DiagDesign[[j]]<-DiagDesign[[j]]*ZZCorrTempTotUpd[,(1+j-currentindex[1])]
                if (j==currentindex[1]) ZZCorrTotUpd<-DiagDesign[[j]]
                else ZZCorrTotUpd<-cbind(ZZCorrTotUpd,DiagDesign[[j]])
                #DiagDesign[[j]]<-model.matrix(~SSC[[j]]-1)
            }
        
            ZZCorrUpd[[i]]<-ZZCorrTotUpd 
        
            if (i==1) ZZCorrDesign<-ZZCorrUpd[[1]]
            else ZZCorrDesign<-cbind(ZZCorrDesign,ZZCorrUpd[[i]])
        }
        q<-sapply(DiagDesign,ncol)
    }
                    
          
        # DiagDesign contains individual design matrix for each random effects - which are now independent !!!!!!!!!!!!! #
        # The models must be fitted jointly if there are correlations between random effects between different models #
    
        # !!!!! It works !!!!! #
        
        # Now we need to add independent design matricies and then create the vector of random effect corresponding to the design matrices #
        # Independent random effects are in ZZIndep indexed by 
        
    # Handle independent random effects here #
    
    if (exists("ZZIndep")) nrandind<-length(ZZIndep)
    else nrandind<-0
    if (exists("ZZCorr")) nrandcor<-length(ZZCorr)
    else nrandcor<-0
      
    Beta<-unlist(BetaList)
    if (nrandind>0) {
        ZZIndepTemp<-list(0)
        for (i in 1:nrandind){
            if (indepModel[i]==length(ModelsDims)) rowsdown<-0
            else rowsdown<-sum(ModelsDims[(indepModel[i]+1):length(ModelsDims)][!is.na(ModelsDims[(indepModel[i]+1):length(ModelsDims)])])
            if (indepModel[i]==1) rowsup<-0
            else rowsup<-sum(ModelsDims[1:(indepModel[i]-1)])
            
            # First expand up if possible #
            ZZIndepTemp[[i]]<-FillZeros(ZZIndep[[i]],nr=(nrow(ZZIndep[[i]])+rowsup),nc=ncol(ZZIndep[[i]]),up=TRUE)
            ZZIndepTemp[[i]]<-FillZeros(ZZIndepTemp[[i]],nr=(nrow(ZZIndepTemp[[i]])+rowsdown),nc=ncol(ZZIndepTemp[[i]]))
            if (!exists("ZZCorrDesign") & i==1) ZZDesign<-ZZIndepTemp[[i]]
            if (exists("ZZCorrDesign") & i==1) ZZDesign<-cbind(ZZCorrDesign,ZZIndepTemp[[1]])
            if (i>1) ZZDesign<-cbind(ZZDesign,ZZIndepTemp[[i]])
            if (i==1) ZZIndepDesign<-ZZIndepTemp[[1]]
            if (i>1) ZZIndepDesign<-cbind(ZZIndepDesign,ZZIndepTemp[[i]])
        }
        
    }
    else if(nrandcor>0) ZZDesign<-ZZCorrDesign
    
    if (nrandind>0) {
        if (!is.null(CorrMat)) q<-c(q,sapply(ZZIndepTemp,ncol))
        else q<-sapply(ZZIndepTemp,ncol)
    }
    
    # Create ZOriginal and TTOriginal #
    # ZOriginal is not well defined here !!!!!! - it must be corrected #
    ZOriginal<-NULL
    if (nrandcor>0) {
        for (i in 1:nrandcor) {# These also dont need to be from different models #
            if (corrModel[i]==length(ModelsDims)) rowsdown<-0
            else rowsdown<-sum(ModelsDims[(corrModel[i]+1):length(ModelsDims)][!is.na(ModelsDims[(corrModel[i]+1):length(ModelsDims)])])
            if (corrModel[i]==1) rowsup<-0
            else rowsup<-sum(ModelsDims[1:(corrModel[i]-1)])
            
            ZZCorrTemp<-FillZeros(ZZCorr[[i]],nr=(nrow(ZZCorr[[i]])+rowsup),nc=ncol(ZZCorr[[i]]),up=TRUE)
            ZZCorrTemp<-FillZeros(ZZCorrTemp,nr=(nrow(ZZCorrTemp)+rowsdown),nc=ncol(ZZCorrTemp))
            if (i==1) ZOriginal<-ZZCorrTemp
            else {
                ZOriginal<-cbind(ZOriginal,ZZCorrTemp)
            }
        }
    }
    
    if (nrandind>0) {
            if (i==1 & is.null(ZOriginal)) ZOriginal<-ZZIndepDesign
            else ZOriginal<-cbind(ZOriginal,ZZIndepDesign)
    }
    
    TTOriginal1<-cbind(X,ZOriginal)
    TTOriginal2<-cbind(matrix(0,ncol(ZOriginal),ptot),diag(ncol(ZOriginal)))
    TTOriginal<-rbind(TTOriginal1,TTOriginal2)
    
    ###################    
    qcum<-cumsum(c(0,q))
    qtot<-qcum[length(qcum)]
    if (is.null(Vstart)) Vstart<-rep(0,sum(q))
          
    
    # index of correlated random effects 
    V<-list(0)
    U<-list(0)
##    print(nrow(ZZDesign))
##    print(ncol(ZZDesign))
##    print(nrow(X))
    TT<-cbind(X,ZZDesign)
    PsiM<-rep(0,sum(q))
    
    if (nrandcor>0){
        for (i in 1:nrandcor){
            TT<-rbind(TT,cbind(matrix(0,q[i],ptot+qcum[i]),diag(q[i]),matrix(0,q[i],qcum[nrandcor+nrandind+1]-qcum[i+1])))
            V[[i]]<-as.matrix(Vstart[(qcum[i]+1):(qcum[i+1])])
            if (i==1) VT<-V[[1]]
            else VT<-c(VT,list(V[[i]]))
            if (RandDistCorr[i]=="Normal") PsiM[(qcum[i]+1):qcum[i+1]]<-0
            if (RandDistCorr[i]=="Gamma")  PsiM[(qcum[i]+1):qcum[i+1]]<-1
            if (RandDistCorr[i]=="IGamma") PsiM[(qcum[i]+1):qcum[i+1]]<-1
            if (RandDistCorr[i]=="Beta")   PsiM[(qcum[i]+1):qcum[i+1]]<-0.5
        }
    }
    
    if (nrandind>0){
        for (i in 1:nrandind){
            TT<-rbind(TT,cbind(matrix(0,q[i+nrandcor],ptot+qcum[i+nrandcor]),diag(q[i+nrandcor]),matrix(0,q[i+nrandcor],qcum[nrandind+nrandcor+1]-qcum[nrandcor+i+1])))
            V[[i+nrandcor]]<-as.matrix(Vstart[(qcum[i+nrandcor]+1):(qcum[i+nrandcor+1])])
            if ((i+nrandcor)==1) VT<-V[[i]]
            else VT<-c(VT,list(V[[i+nrandcor]]))
            if (RandDistIndep[i]=="Normal") PsiM[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]<-0
            if (RandDistIndep[i]=="Gamma")  PsiM[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]<-1
            if (RandDistIndep[i]=="IGamma") PsiM[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]<-1
            if (RandDistIndep[i]=="Beta")   PsiM[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]<-0.5
        }
    }  
    
    # OFFSET and Binomial denominators management - i am not sure how to do it - keep it for later while we develop the code # 
    # The offset needs to be done for the whole vector corresponding to Y - such as all models are considered Poisson or Binomial -begrijp je? #
    Z<-TT[(1:sum(ModelsDims)),((ptot+1):ncol(TT))]

    B<-0
    if (is.null(OFFSETList)) B<-rep(1,nrow(Y))
    else {
        for (i in 1:nModels){
            if (length(OFFSETList[[i]])==1 & OFFSETList[[i]]==0) B[(cModelsDims[i]+1):(cModelsDims[i+1])]<-1
            else B[(cModelsDims[i]+1):(cModelsDims[i+1])]<-OFFSETList[[i]]       
        }
    }
    
    
    if (nrandcor>0) DRgamma<-rep(0,ncol(DDRCorr))
    else DRgamma<-NULL
    if (nrandind>0) DRgamma<-c(DRgamma,DRgammaIndep)
        
        
    Iteration<-0
    Convergence<-100

     while (Convergence>CONV){ 
#     for (iii in 1:1) {   
        Iteration<-Iteration+1
#        if (Info) cat("\n Iteration: ",Iteration,"     Convergence: ",Convergence,"\n")
        
        ###############################################
        # PLAN:                                       #
        # Update Mean Structure                       #
        # Laplace approximation to the mean structure #
        # Variances of random components              #
        # Overdispersions                             #
        # Correlations                                #
        ###############################################
        
        MeanParmsLast<-c(Beta,unlist(VT))
        # Lambda of correlated random effects is equal to one - these after cholesky transformation are denoted LambdaChol#
        if (nrandcor>0) {
            DRgamma[1:cumqcorr[length(cumqcorr)]]<-0
        }
        Lambda<-exp(DDR%*%DRgamma)

        # Overdispersion - now for all models overdispersion is coded in one matrix DDY and one vector DYgamma FOR ALL MODELS AT THE SAME TIME #
        Phi<-exp(DDY%*%DYgamma)
        GammaMvec<-c(Phi,Lambda)
        
        # Mean values #
        eta<-TT[1:ntot,]%*%as.matrix(c(Beta,unlist(VT)))

        # Here different link is applied to different model #
        mu<-0
        Wvec<-0
        dWdmu<-0
        dmudeta<-0
        d2Wdmu2<-0
        for (i in 1:nModels){
            mu[(cModelsDims[i]+1):cModelsDims[i+1]]<-B[(cModelsDims[i]+1):cModelsDims[i+1]]*InvLinkY(eta[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]])
            Wvec[(cModelsDims[i]+1):cModelsDims[i+1]]<-Wmatgen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
            dWdmu[(cModelsDims[i]+1):cModelsDims[i+1]]<-dWdmugen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
            d2Wdmu2[(cModelsDims[i]+1):cModelsDims[i+1]]<-d2Wdmu2gen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
            dmudeta[(cModelsDims[i]+1):cModelsDims[i+1]]<-dmudetagen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        }
     
        # So far the Link and Dist are equivalent as only canonical links are allowed #
        WR<-list(0)
        UT<-0
        WRT<-list(0)
        if (nrandcor+nrandind>0) {
            for (i in 1:(nrandcor+nrandind)) {
                U[[i]]<-LinkR(V[[i]],RandDist[i])
                if (i==1) UT<-U[[1]]
                else UT<-c(UT,U[[i]])
                WR[[i]]<-WRVC(U[[i]],RandDist[i])
                if (i==1) WRT<-WR[[1]]
                else WRT<-c(WRT,WR[[i]])
            }
        }        
        
        WTotvec<-c(Wvec,WRT)
                
        ISIGMAMvec<-as.vector((1/GammaMvec)*WTotvec)       
       
        #dAdmu<-dAdmugen(mu,Link)
        #Amat<-Amatgen(mu,Link)
        # Adjustment computation for the Laplace Approximation to the mean #
        
        # Now a Laplace Approximation to the mean can be used for one model but not for another therefore LaplaceFixed is a vector #
        
        # We will not make a LAPFIX by model, instead we adjust all at the same time:
        #   Those models which are LAPFIX=FALSE need to have Design matrix in this part Z equal to zero for random effects corresponding to these models 
        #   The adjustment terms are of length n1+n2+.... so for all the models, if one model is LAPFIX=FALSE then the adjustemt terms for those models are evenly 
        #   Redistributed among the rest of the models which are LAPFIX=TRUE
        
        Scorr<-rep(0,sum(ModelsDims))
        
        ################################################
        # Consider the revision of the procedure below #
        # By adding the factor random effects          #
        ################################################
            
        # There are two options : one is that there are too many adjustments than the ys - so they have to be redistributed
        #                       : two is that rows are deleted and there are as many adjustments as the ys but instead all the dimension change
        #                       : total matrices must be used in computation of the derivitive vhatbeta
             
        # PART OF CREATING DESIGN MATRICES FOR LAPFIX CAN BE MOVED OUT OF THE ITERATIONS AS THE SAME MATRICES WILL BE USED - BUT ONLY CORRELATION CHANGES THE 
        # CORRELATED PART - THINK HOW TO DO IT      
        
        CorrTerms<-0
        if (any(LaplaceFixed)==TRUE){
            
            # Exlcude from the Z matrix the collumns which are not integreatd out - as this model is not LAPFIX=TRUE #
            # We need an index to denote which random effects are integrated out #
            ZLapFix<-Z
            VLapFix<-unlist(VT)
            qLapFix<-q
            ISIGMAMvecLapFix<-ISIGMAMvec
            WvecLapFix<-Wvec
            dWdmuLapFix<-dWdmu
            dmudetaLapFix<-dmudeta
            integratedModels<-seq(1,nModels)
            ModelsDimsLapFix<-ModelsDims
            
            if (nrandcor>0) {
                for (i in 1:nrandcor){ 
                    if (LaplaceFixed[corrModel[i]]==FALSE) {
                        ZLapFix[,(qcum[i]+1):(qcum[i+1])]<-0
                        VLapFix[(qcum[i]+1):(qcum[i+1])]<-NA
                        qLapFix[i]<-NA
                        ISIGMAMvecLapFix[(ntot+qcum[i]+1):(ntot+qcum[i+1])]<-NA
                    }
                }
            }
            
            if (nrandind>0) {
                for (i in 1:nrandind){
                    if (LaplaceFixed[indepModel[i]]==FALSE) {
                        ZLapFix[,(qcum[i+nrandcor]+1):(qcum[i+nrandcor+1])]<-0
                        VLapFix[(qcum[i+nrandcor]+1):(qcum[i+nrandcor+1])]<-NA
                        qLapFix[i+nrandcor]<-NA
                        ISIGMAMvecLapFix[(ntot+qcum[i+nrandcor]+1):(ntot+qcum[i+nrandcor+1])]<-NA
                    }
                }
            }
            
            nrandint<-0
            integratedIndex<-NA
            if (nrandcor>0) {
                for (i in 1:nrandcor){
                    if (LaplaceFixed[corrModel[i]]==TRUE) {
                        nrandint<-nrandint+1
                        if (is.na(integratedIndex[1])) integratedIndex<-i
                        else integratedIndex<-c(integratedIndex,i)
                    }
                }
            }
            if (nrandind>0) {
                for (i in 1:nrandind){
                    if (LaplaceFixed[indepModel[i]]==TRUE) {
                        nrandint<-nrandint+1
                        if (is.na(integratedIndex[1])) integratedIndex<-i
                        else integratedIndex<-c(integratedIndex,i+nrandcor) # integratedIndex contains the number of random effect which is integrated out its position in original design matrix#
                    }
                }
            }
            # Remove all the columns which are of zero - but question should we keep track which columns are removed ? #
            
            ZLapFix<-ZLapFix[,(apply(ZLapFix,2,sum)!=0)]
            VLapFix<-VLapFix[!is.na(VLapFix)]
            qLapFix<-qLapFix[!is.na(qLapFix)]
            qcumLapFix<-cumsum(c(0,qLapFix))
            
            
            # Remove Rows if are all zeros - and create an index which rows are removed - note that we remove the rows for Y response together with #
            rowLapFix<-rep(1,nModels) # This index will determine which rows are deleted for which model if index is set to zero it means that the rows for that model 
                                        # are deleted
            obsLapFix<-rep(1,ntot)   # This vector will select which adjustments are assigned to ys and which need to be redistributed at the moment i am doing it, it is
                                        # still unclear how it will work
            
            for (i in 1:nModels) {
                if (all(apply(ZLapFix[(cModelsDims[i]+1):cModelsDims[i+1],],1,sum)==0)) {
                    integratedModels[i]<-NA
                    rowLapFix[i]<-0
                    ISIGMAMvecLapFix[(cModelsDims[i]+1):cModelsDims[i+1]]<-NA # This removes the main effects - we still need to remove the random effects matrix #
                    WvecLapFix[(cModelsDims[i]+1):cModelsDims[i+1]]<-NA
                    dWdmuLapFix[(cModelsDims[i]+1):cModelsDims[i+1]]<-NA
                    dmudetaLapFix[(cModelsDims[i]+1):cModelsDims[i+1]]<-NA
                    ModelsDimsLapFix[i]<-NA
                    
                    if (nrandcor>0){
                        for (j in 1:nrandcor){
                            if (corrModel[j]==i) ISIGMAMvecLapFix[(ntot+qcum[j]+1):(ntot+qcum[j+1])]<-NA
                        }
                    }
                    if (nrandind>0){
                        for (j in 1:nrandind){
                            if (indepModel[j]==i) ISIGMAMvecLapFix[(ntot+qcum[j+nrandcor]+1):(ntot+qcum[j+nrandcor+1])]<-NA
                        }
                    }
                }
                if (LaplaceFixed[i]==FALSE) obsLapFix[(cModelsDims[i]+1):cModelsDims[i+1]]<-0
                
            }
            
            ZLapFix<-ZLapFix[apply(ZLapFix,1,sum)!=0,]
            ISIGMAMvecLapFix<-ISIGMAMvecLapFix[!is.na(ISIGMAMvecLapFix)]
            WvecLapFix<-WvecLapFix[!is.na(WvecLapFix)] 
            dWdmuLapFix<-dWdmuLapFix[!is.na(dWdmuLapFix)]
            dmudetaLapFix<-dmudetaLapFix[!is.na(dmudetaLapFix)]
            ntotLapFix<-length(WvecLapFix)
            integratedModels<-integratedModels[!is.na(integratedModels)]
            ModelsDimsLapFix<-ModelsDimsLapFix[!is.na(ModelsDimsLapFix)]
            cModelsDimsLapFix<-cumsum(c(0,ModelsDimsLapFix))
            
            nintMod<-length(integratedModels)
            # We need to separate the design matrix for the derivative with respect to beta with the design which is used in the trace part #
            
            TT2<-TT[,(ptot+1):(ptot+sum(q))]
            TT2LapFix<-rbind(ZLapFix,diag(ncol(ZLapFix)))            
            
            # The function Solver Short needs to be reprogrammed - but how can we invert a matrix if there is a zero row or zero collumn #
            # One thing is to use a generalized inverse, however other thing is to - do we need zeros in the inverse function #
            
            # This is for the derivative #
            OUT1<-SolverShort(ISIGMAMvec,Z) 
            OUTLapFix<-SolverShort(ISIGMAMvecLapFix,ZLapFix) # This needs to be adjusted to accomodate the zero rows, than the diagonal matrix is reduced #
            # A new argument Z is added.as the design matrix in this part might be different than the general Z
            #print("DimTT2");print(dim(TT2));print("Dim2");print(length(rep(ISIGMAMvec,each=nrow(t(TT2)))))
            INV1<-OUT1$Inverse              # These are used for derivative of random effects #
            DiagPP2<-OUT1$DiagPP2
            INV1LapFix<-OUTLapFix$Inverse   # These are used for the correction factor determinant over integrated random effects #
            DiagPP2LapFix<-OUTLapFix$DiagPP2 # This will not work because ntot and dimensions are given globally we need to change it !#
            
            rm(OUT1)    
            rm(OUTLapFix)            

            MOD<-INV1%*%(t(Z)*rep((1/Phi),each=ncol(Z))) # This is for random effects derivative
            ADJDER1<-list(0)
            ADJDER2<-list(0)
            
            # How to modify this? #
            # Now we need to iterate over the random effects which are integrated out #
            for (i in 1:nrandint){
                ADJDER1[[i]]<--ZLapFix[,(qcumLapFix[i]+1):qcumLapFix[i+1]]%*%(MOD[(qcum[integratedIndex[i]]+1):(qcum[integratedIndex[i]+1]),])#*rep(Wvec,each=nrow(MOD[(qcum[i]+1):(qcum[i+1]),]))) 
                ADJDER2[[i]]<--MOD[(qcum[integratedIndex[i]]+1):(qcum[integratedIndex[i]+1]),]
            }        

            # Computing correction quantities for the Laplace Approximation of fixed effects #
            # This here gets difficult ntot and ntotLapFix - how to make sure we do the correct thing #
            CorrTermsTemp<-list(0)
            CorrTermsTemp[[1]]<-as.vector(DiagPP2LapFix[1:ntotLapFix])*as.vector(1/WvecLapFix)*as.vector(1/WvecLapFix)*(dWdmuLapFix)*dmudetaLapFix
            # We need to scale CorrTerms[[1]] onto the whole ntot #
            CorrTerms<-list(0)
            CorrTerms[[1]]<-rep(0,ntot)
            for (i in 1:nintMod){
                CorrTerms[[1]][(cModelsDims[integratedModels[i]]+1):cModelsDims[integratedModels[i]+1]]<-CorrTermsTemp[[1]][(cModelsDimsLapFix[i]+1):cModelsDimsLapFix[i+1]]
            }
            
            for (i in 1:nrandint){
                ADJ1<-rep(0,ntot)
                ADJ2<-rep(0,ntot)
                ADJ1<-t(ADJDER1[[i]])%*%(as.vector(DiagPP2LapFix[1:ntotLapFix])*as.vector(1/WvecLapFix)*as.vector(dWdmuLapFix)*as.vector(dmudetaLapFix)) # Check this one
                if (RandDist[integratedIndex[i]]=="Gamma") ADJ2<-t(ADJDER2[[i]])%*%as.vector(DiagPP2LapFix[(ntotLapFix+qcumLapFix[i]+1):(ntotLapFix+qcumLapFix[1+i])])
                if (RandDist[integratedIndex[i]]=="IGamma") ADJ2<-t(ADJDER2[[i]])%*%(as.vector(DiagPP2LapFix[(ntotLapFix+qcumLapFix[i]+1):(ntotLapFix+qcumLapFix[1+i])])*as.vector(2*U[[integratedIndex[i]]]))
                if (RandDist[integratedIndex[i]]=="Beta") ADJ2<-t(ADJDER2[[i]])%*%(as.vector(DiagPP2LapFix[(ntotLapFix+qcumLapFix[i]+1):(ntotLapFix+qcumLapFix[1+i])])*as.vector(1-2*U[[integratedIndex[i]]]))
                CorrTerms<-c(CorrTerms,list(ADJ1,ADJ2))
            }
            
            # Here we need to have in mind two scenarions, when there is cholesky correlation therefore design over n1+n2 but integrated only the one model #
            # And scenario two when only one model is integrated but there are no cholesky correlation #
            CorrTermsLength<-length(CorrTerms)
            CorrTerms<-as.matrix(unlist(CorrTerms))
          
            dim(CorrTerms)<-c(ntot,CorrTermsLength)

            Scorr<-0.5*Phi*dmudeta*apply(CorrTerms,1,sum)            
            
            # Now we need to take out the corrections for which ys are not adjusted and split them evenly over those who are #
            
            
        }
    
        Ystar<-Y-Scorr  
        zmain<-eta+(Ystar-mu)/dmudeta

        PsiMstar<-PsiM+(Lambda*crossprod(Z,as.matrix(Wvec*Scorr*(1/Phi)/dmudeta)))
        zrand<-list(0)
        if (nrandcor>0){
            for (i in 1:nrandcor){
                zrand[[i]]<-V[[i]]+(PsiMstar[(qcum[i]+1):qcum[i+1]]-U[[i]])/WR[[i]]
            }
        }
        if (nrandind>0){
            for (i in 1:nrandind){
                zrand[[i+nrandcor]]<-V[[i+nrandcor]]+(PsiMstar[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]-U[[i+nrandcor]])/WR[[i+nrandcor]]
            }
        }
        zrand<-as.matrix(unlist(zrand))
         
        zTot<-as.matrix(c(zmain,zrand))
         
        # Updating Equations #
        MeanParmsLast<-c(Beta,unlist(VT))
        
        OUT1<-SolverLong(ISIGMAMvec,zTot)
        MeanParms<-OUT1$NewParms
        DiagPMAT<-OUT1$DiagPMAT
        rm(OUT1)

        if (exists("INV1")) rm(INV1)
         
         #print("Block 2");print(CPTISIGMAM);#print("EISIGMAMvec");#print(EISIGMAMvec);print("zTot");print(zTot)
         
        Beta<-MeanParms[1:ptot]
        if (nrandcor>0){
            for (i in 1:nrandcor){
                V[[i]]<-MeanParms[(ptot+qcum[i]+1):(ptot+qcum[i+1])]
                if (i==1) VT<-V[[i]]
                else VT<-c(VT,V[[i]])
            }
        }
         
        if (nrandind>0){
            for (i in 1:nrandind){
                V[[i+nrandcor]]<-MeanParms[(ptot+qcum[i+nrandcor]+1):(ptot+qcum[i+nrandcor+1])]
                if (i==1) VT<-V[[i+nrandcor]]
                else VT<-c(VT,V[[i+nrandcor]])
            }
        }
        
        Convergence<-sum(abs(MeanParms-MeanParmsLast)) 
##        if (DEBUG==TRUE) print("Convergence Mean");print(Convergence)
        # Now we move to the estimation of dispersion and overdispersion given the correlation parameters rho based on transformed multivariat random effects to the 
        # independence scenario - so we continue the procedure, first dispersion for each random effect, then overdispersion for each model 
        # After that we reestimate the correlation and construct updated variance covariance matrix this ends the joint model
        
       ###############################
       ##### Variance Components #####
       ###############################
       if (EstimateVariances==TRUE) {
       # Reevaluation of mean and u #
       eta<-TT[1:ntot,]%*%as.matrix(c(Beta,unlist(VT)))
       # Here different link is applied to different model #
       for (i in 1:nModels){
            mu[(cModelsDims[i]+1):cModelsDims[i+1]]<-B[(cModelsDims[i]+1):cModelsDims[i+1]]*InvLinkY(eta[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]])
            Wvec[(cModelsDims[i]+1):cModelsDims[i+1]]<-Wmatgen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
            dWdmu[(cModelsDims[i]+1):cModelsDims[i+1]]<-dWdmugen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
            dmudeta[(cModelsDims[i]+1):cModelsDims[i+1]]<-dmudetagen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        }
        if (DEBUG==TRUE) {
            muout <- mu
            Wvecout <- Wvec
            dWdmuout <- dWdmu
            dmudetaout <- dmudeta 

        }
        # So far the Link and Dist are equivalent as only canonical links are allowed #
        WR<-list(0)
        UT<-0
        WRT<-0
        dWRdu<-0
        if (nrandcor+nrandind>0) {
            for (i in 1:(nrandcor+nrandind)) {
                U[[i]]<-LinkR(V[[i]],RandDist[i])
                if (i==1) UT<-U[[1]]
                else UT<-c(UT,U[[i]])
                WR[[i]]<-WRVC(U[[i]],RandDist[i])
                if (i==1) WRT<-WR[[1]]
                else WRT<-c(WRT,WR[[i]])
                if (i==1) dWRdu<-DWRDU(U[[i]],RandDist[i])
                else dWRdu<-c(dWRdu,DWRDU(U[[i]],RandDist[i]))
            }
        }        
        
        WTotvec<-c(Wvec,WRT)
                
        ISIGMAMvec<-as.vector((1/GammaMvec)*WTotvec)      
        
        #Computing Deviances : this must be adjusted for correlated and independent random effects #
        
        DevianceRand<-list(0)
        if (nrandcor>0) {
            for (i in 1:nrandcor){
                DevianceRand[[i]]<-DevRand(V[[i]],U[[i]],RandDist[i])
            }
        }
        if (nrandind>0) { 
            for (i in 1:nrandind){
                DevianceRand[[i+nrandcor]]<-DevRand(V[[i+nrandcor]],U[[i+nrandcor]],RandDist[i+nrandcor])
            }
        }
                
        
        # Truncated Computations #
        #MTheta<-1-exp(-mu)
        #M1Theta<-exp(-mu)*mu
        #M2Theta<-exp(-mu)*mu*(1-mu)
        #M3Theta<-M2Theta*(1-mu)-mu*M1Theta
        #WTildevec<-as.vector(Wvec+((M2Theta/MTheta)-(M1Theta/MTheta)^2))

        DVhatDlambda<-list(0)
          
        # ERROR - This part needs to be adjusted !!!!!!! #
        
        # OLD CODE FOR DERIVATIVE #
        if (nrandcor>0) {          
            for (i in 1:nrandcor){
                DVhatDlambda[[i]]<--solve((t(Z[,(qcum[i]+1):qcum[i+1]])*rep(Wvec/Phi,each=ncol(Z[,(qcum[i]+1):qcum[i+1]])))%*%Z[,(qcum[i]+1):qcum[i+1]]+diag(WR[[i]]/Lambda[(1+qcum[i]):qcum[i+1]]))%*%(PsiM[(qcum[i]+1):(qcum[i+1])]-U[[i]])/(Lambda[(1+qcum[i]):qcum[i+1]]^2)
            if (i==1) derold<-DVhatDlambda[[i]]
            }
        }
        
        if (nrandind>0) {
            for (i in 1:nrandind){
                DVhatDlambda[[i+nrandcor]]<--solve((t(Z[,(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]])*rep(Wvec/Phi,each=ncol(Z[,(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]])))%*%Z[,(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]+diag(WR[[i+nrandcor]]/Lambda[(1+qcum[i+nrandcor]):qcum[i+nrandcor+1]]))%*%(PsiM[(qcum[i+nrandcor]+1):(qcum[i+nrandcor+1])]-U[[i+nrandcor]])/(Lambda[(1+qcum[i+nrandcor]):qcum[i+nrandcor+1]]^2)
            }
        }
        
        # New code for derivative #
        # We need to create the lambda vector #
        #if (nrandcor>0) {
        #    for (i in 1:nrandcor) {
        #        LambdaTemp<-rep(0,length(Lambda))
        #        LambdaTemp[(qcum[i]+1):qcum[i+1]]<-1
        #        DVhatDlambda[[i]]<--solve(t(Z*as.vector(Wvec/Phi))%*%Z+diag(as.vector(WRT/Lambda)))%*%(((PsiM-UT)/(Lambda^2))*LambdaTemp)
        #        if (i==1) dernew<-DVhatDlambda[[i]]
        #    }
        #}
        
        #if (nrandind>0) {
        #    for (i in 1:nrandind) {
        #        LambdaTemp<-rep(0,length(Lambda))
        #        LambdaTemp[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]<-Lambda[(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]
        #        DVhatDlambda[[i+nrandcor]]<--solve(t(Z*as.vector(Wvec/Phi))%*%Z+diag(as.vector(WRT/Lambda)))%*%(((PsiM-UT)/(Lambda^2))*LambdaTemp)
        #    }
        #}        
        
            
        #DWTildeDthetavec<-as.vector(WTildevec-(M2Theta/MTheta)+(M3Theta/MTheta)-3*((M1Theta*M2Theta)/(MTheta^2))+((M1Theta^2)/(MTheta^2))+2*((M1Theta^3)/(MTheta^3)))
        #DWTildeDmuvec<-as.vector(DWTildeDthetavec/mu)
              
        qmod<-list(0)
        qCur<-list(0)
        qrr<-list(0)
        if (nrandcor>0) {
            for (i in 1:nrandcor){
                SSCur<-rep(0,ntot+sum(q))
                SSCur[1:ntot]<-SSC[[i]]
                SSCur[(ntot+1+qcum[i]):(ntot+qcum[i+1])]<-(1:(qcum[i+1]-qcum[i]))
            
                DlogWDloglambda<-matrix(0,ntot+sum(q),ntot+sum(q))
                ##### This here can be replaced the diagonal matrix dlogwdloglambda #####
                DlogWDloglambda[1:ntot,1:ntot]<-diag(as.vector(1/Wvec)*as.vector(dWdmu)*as.vector(dmudeta)*(as.vector(Z[,(1+qcum[i]):qcum[i+1]]%*%(DVhatDlambda[[i]]*Lambda[(1+qcum[i]):qcum[i+1]])))) # RRRRR
                DlogWDloglambda[(ntot+1+qcum[i]):(ntot+qcum[i+1]),(ntot+1+qcum[i]):(ntot+qcum[i+1])]<-diag(DWRDU(U[[i]],RandDist[i])*(as.vector(DVhatDlambda[[i]])*as.vector(Lambda[(1+qcum[i]):qcum[i+1]]))) # RRRRR
                qmod[[i]]<-DiagPMAT*diag(DlogWDloglambda)
                qCur[[i]]<-cbind(qmod[[i]],SSCur)

                qCur[[i]]<-tapply(qCur[[i]][,1],qCur[[i]][,2],sum)

                qCur[[i]]<-qCur[[i]][row.names(qCur[[i]])!="0"]

                qrr[[i]]<-DiagPMAT[(ntot+1+qcum[i]):(ntot+qcum[i+1])]
##                if (DEBUG) {print("Basis qrr");print(qrr);}
                qrr[[i]]<-qrr[[i]]-qCur[[i]]
                # Correction to estimate the true likelihood instead of EQL #
                if (RandDist[i]=="Gamma") qrr[[i]]<-qrr[[i]]+1+2*(log(Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])+2*(digamma(1/Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])
                if (RandDist[i]=="IGamma") qrr[[i]]<-qrr[[i]]+1+(2/Lambda[(1+qcum[i]):(qcum[i+1])])-2*(log(1/Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])-2*((1+Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])+2*(digamma(1+(1/Lambda[(1+qcum[i]):(qcum[i+1])]))/Lambda[(1+qcum[i]):(qcum[i+1])])
                if (RandDist[i]=="Beta")  qrr[[i]]<-qrr[[i]]+1-2*(digamma(1/Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])+2*(digamma(1/(2*Lambda[(1+qcum[i]):(qcum[i+1])]))/Lambda[(1+qcum[i]):(qcum[i+1])])+log(4)/Lambda[(1+qcum[i]):(qcum[i+1])]

                # Applying the correction for the deviances #
                DevianceRand[[i]]<-DevianceRand[[i]]/(1-qrr[[i]]) 
            
            }
        }
        
       if (nrandind>0) {
            for (i in 1:nrandind){
                SSCur<-rep(0,ntot+sum(q))
                SSCur[1:ntot]<-SSIndep[[i]]
                SSCur[(ntot+1+qcum[i+nrandcor]):(ntot+qcum[i+nrandcor+1])]<-(1:(qcum[i+nrandcor+1]-qcum[i+nrandcor]))
            
                DlogWDloglambda<-matrix(0,ntot+sum(q),ntot+sum(q))
                ##### This here can be replaced the diagonal matrix dlogwdloglambda #####
                DlogWDloglambda[1:ntot,1:ntot]<-diag(as.vector(1/Wvec)*as.vector(dWdmu)*as.vector(dmudeta)*(as.vector(Z[,(qcum[i+nrandcor]+1):qcum[i+nrandcor+1]]%*%(DVhatDlambda[[i+nrandcor]]*Lambda[(1+qcum[i+nrandcor]):qcum[i+nrandcor+1]])))) # RRRRR
                DlogWDloglambda[(ntot+1+qcum[i+nrandcor]):(ntot+qcum[i+nrandcor+1]),(ntot+1+qcum[i+nrandcor]):(ntot+qcum[i+nrandcor+1])]<-diag(DWRDU(U[[i+nrandcor]],RandDist[i+nrandcor])*(as.vector(DVhatDlambda[[i+nrandcor]])*as.vector(Lambda[(1+qcum[i+nrandcor]):qcum[i+nrandcor+1]]))) # RRRRR
                qmod[[i+nrandcor]]<-DiagPMAT*diag(DlogWDloglambda)
                qCur[[i+nrandcor]]<-cbind(qmod[[i+nrandcor]],SSCur)

                qCur[[i+nrandcor]]<-tapply(qCur[[i+nrandcor]][,1],qCur[[i+nrandcor]][,2],sum)

                qCur[[i+nrandcor]]<-qCur[[i+nrandcor]][row.names(qCur[[i+nrandcor]])!="0"]

                qrr[[i+nrandcor]]<-DiagPMAT[(ntot+1+qcum[i+nrandcor]):(ntot+qcum[i+nrandcor+1])]
##                if (DEBUG) {print("Basis qrr");print(qrr);}
                qrr[[i+nrandcor]]<-qrr[[i+nrandcor]]-qCur[[i+nrandcor]]
                # Correction to estimate the true likelihood instead of EQL #
                if (RandDist[i+nrandcor]=="Gamma") qrr[[i+nrandcor]]<-qrr[[i+nrandcor]]+1+2*(log(Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])+2*(digamma(1/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])
                if (RandDist[i+nrandcor]=="IGamma") qrr[[i+nrandcor]]<-qrr[[i+nrandcor]]+1+(2/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])-2*(log(1/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])-2*((1+Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])+2*(digamma(1+(1/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])]))/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])
                if (RandDist[i+nrandcor]=="Beta")  qrr[[i+nrandcor]]<-qrr[[i+nrandcor]]+1-2*(digamma(1/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])+2*(digamma(1/(2*Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])]))/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])])+log(4)/Lambda[(1+qcum[i+nrandcor]):(qcum[i+nrandcor+1])]

                # Applying the correction for the deviances #
                DevianceRand[[i+nrandcor]]<-DevianceRand[[i+nrandcor]]/(1-qrr[[i+nrandcor]]) 
            
            }
        }
        
##        if (DEBUG) {print("Qrr");print(qrr);}
        ######################################
        ##### Gamma model for dispersion #####
        ######################################

        # To speed up we can separate the models for dispersion by component instead of inverting the whole matrix #
        # However if there are not so many dispersion components it should not be a problem #
            
        invSigmaGammaR<-list(0)
        if (nrandcor>0){
            for (i in 1:nrandcor){ 
                invSigmaGammaR[[i]]<-((1-qrr[[i]])/4)
            }
        }
        
        if (nrandind>0){
            for (i in 1:nrandind){ 
                invSigmaGammaR[[i+nrandcor]]<-((1-qrr[[i+nrandcor]])/4)
            }
        }
                
        muGammaR<-Lambda
        
        oldDRgamma<-DRgamma
        
        ksiR<-DDR%*%DRgamma

        DevianceRand<-unlist(DevianceRand)
        
        ZRresp<-ksiR+(DevianceRand-muGammaR)/muGammaR
        
        invSigmaGammaR<-diag(unlist(invSigmaGammaR))
        if (nrandcor==0) cumqcorr <- 0
        OldIndepgamma <- DRgamma[(cumqcorr[length(cumqcorr)]+1):length(DRgamma)]
        DRgamma<-solve(crossprod(DDR,invSigmaGammaR)%*%DDR,crossprod(DDR,invSigmaGammaR)%*%ZRresp)
        StdErrDRgamma<-sqrt(diag(solve(crossprod(DDR,invSigmaGammaR)%*%DDR)))
        NewIndepgamma <-DRgamma[(cumqcorr[length(cumqcorr)]+1):length(DRgamma)]
        if (DEBUG==TRUE){
##            print("DRgamma");print(DRgamma)
##            print("Beta");print(Beta)
##            print("VS");print(VT)
        }    
        #if (DEBUG) {print("Phi End");print(Phi)}
        if (nrandind>0) Convergence <- Convergence + sum(abs(OldIndepgamma-NewIndepgamma))
        
        
        # Maybe the piece of code below should be moved after the correlation is estimated #
        ####################################################################################################################
        ##### In this part we redesign the Z matrix - update with new correlation and new variances for the correlated #####
        ####################################################################################################################
        
        # We need to add the construction of correlation matrix after the correlations have been upadted in front of this piece code #
        
        if (nrandcor>0){
            for (i in 1:length(CorrMat)){
                currentGamma<-DRgamma[(cumqcorr[i]+1):cumqcorr[i+1]]
                DiagVar<-diag(exp(currentGamma))
##                if (DEBUG) print("Correlation Matrix");print(DiagVar)
                # Apply inverse cholesky #
                currentMat<-CholeskyMatrices[[i]]%*%DiagVar%*%t(CholeskyMatrices[[i]])
                # Now update with correlations #
                currentMat<-diag(currentMat)
                currentMat<-sqrt(currentMat)*t(CorrMatOut[[i]]*sqrt(currentMat))
                CholeskyMatrices[[i]]<-t(chol(currentMat))
                # Now current mat is the the new varcov matrix with the same correlation as before #
                ZZCorrTempTotUpd<-ZZShort[[i]]%*%CholeskyMatrices[[i]]
                currentindex<-seq(cumqcorr[i]+1,cumqcorr[i+1])
                
                for (j in currentindex) {
                    DiagDesign[[j]]<-model.matrix(~SSC[[j]]-1)
                    DiagDesign[[j]]<-DiagDesign[[j]]*ZZCorrTempTotUpd[,(1+j-currentindex[1])]
                    if (j==currentindex[1]) ZZCorrTotUpd<-DiagDesign[[j]]
                    else ZZCorrTotUpd<-cbind(ZZCorrTotUpd,DiagDesign[[j]])
                }
                
                ZZCorrUpd[[i]]<-ZZCorrTotUpd
                
                if (i==1) ZZCorrDesign<-ZZCorrUpd[[1]]
                else ZZCorrDesign<-cbind(ZZCorrDesign,ZZCorrUpd[[i]])        
            }
        }
         
        if (nrandind>0) {
            ZZIndepTemp<-list(0)
            for (i in 1:nrandind){
                if (indepModel[i]==length(ModelsDims)) rowsdown<-0
                else rowsdown<-sum(ModelsDims[(indepModel[i]+1):length(ModelsDims)])
                if (indepModel[i]==1) rowsup<-0
                else rowsup<-sum(ModelsDims[1:(indepModel[i]-1)])
            
                # First expand up if possible #
                ZZIndepTemp[[i]]<-FillZeros(ZZIndep[[i]],nr=(nrow(ZZIndep[[i]])+rowsup),nc=ncol(ZZIndep[[i]]),up=TRUE)
                ZZIndepTemp[[i]]<-FillZeros(ZZIndepTemp[[i]],nr=(nrow(ZZIndepTemp[[i]])+rowsdown),nc=ncol(ZZIndepTemp[[i]]))
                if (!exists("ZZCorrDesign") & i==1) ZZDesign<-ZZIndepTemp[[i]]
                if (exists("ZZCorrDesign") & i==1) ZZDesign<-cbind(ZZCorrDesign,ZZIndepTemp[[1]])
                if (i>1) ZZDesign<-cbind(ZZDesign,ZZIndepTemp[[i]])
            }
        
        }
        else if(nrandcor>0) ZZDesign<-ZZCorrDesign
        TT<-cbind(X,ZZDesign)
        if (nrandcor>0){
            for (i in 1:nrandcor){
                TT<-rbind(TT,cbind(matrix(0,q[i],ptot+qcum[i]),diag(q[i]),matrix(0,q[i],qcum[nrandcor+nrandind+1]-qcum[i+1])))
            }
        }
    
        if (nrandind>0){
            for (i in 1:nrandind){
                TT<-rbind(TT,cbind(matrix(0,q[i+nrandcor],ptot+qcum[i+nrandcor]),diag(q[i+nrandcor]),matrix(0,q[i+nrandcor],qcum[nrandind+nrandcor+1]-qcum[nrandcor+i+1])))
            }
        }  
    
        Z<-TT[(1:sum(ModelsDims)),((ptot+1):ncol(TT))]         
        
        # Koniec respecyfikacji macierzy efektow losowych w czesci odpowiadajacej efektom skorelowanym #
        #Convergence<-Convergence+sum(abs(DRgamma-oldDRgamma))     
        }
        ##########################
        ##### Overdispersion #####
        ##########################
        
        # now work on overdispersion #
        # we need to select the parameters which are supposed to be EstimateOverDisp=TRUE not all models have estimated overdispersion true #
        
        if (any(EstimateOverDisp)){
            
            # First select the matrix which will be used for the estimation #
            for (i in 1:nModels){
                if (EstimateOverDisp[i]==TRUE) {
                    if (i==1) indexODEstim<-seq(cModelsDims[i]+1,cModelsDims[i+1])
                    else indexODEstim<-c(indexODEstim,seq(cModelsDims[i]+1,cModelsDims[i+1]))
                    ntotODEstim<-length(indexODEstim)                       
                }
            }
                              
            #PMAT<-TT%*%solve(CPTISIGMAM%*%TT)%*%CPTISIGMAM     # This stays the same - just a selection need to be done of the diagonals corresponding to the models which are #
                                                              # EstimOverDisp = TRUE #
              
            qrrO<-rep(0,ntotODEstim)
                   
            # Applying the correction for the deviances #
            DevianceResp<-rep(0,ntotODEstim)
            DiagPMATODEstim<-DiagPMAT
            for (i in 1:nModels){
                if (EstimateOverDisp[i]==TRUE){
                    DevianceRespTemp<-rep(0,cModelsDims[i+1]-cModelsDims[i])
                    YTemp<-YList[[i]]
                    BTemp<-B[(cModelsDims[i]+1):cModelsDims[i+1]]
                    muTemp<-mu[(cModelsDims[i]+1):cModelsDims[i+1]]
                    PhiTemp<-Phi[(cModelsDims[i]+1):cModelsDims[i+1]]
                    if (RespDist[i]=="Binomial") {
                        DevianceRespTemp[YTemp!=0 & YTemp!=BTemp]<-2*(YTemp[YTemp!=0 & YTemp!=BTemp]*log(YTemp[YTemp!=0 & YTemp!=BTemp]/muTemp[YTemp!=0 & YTemp!=BTemp])-(YTemp[YTemp!=0 & YTemp!=BTemp]-BTemp[YTemp!=0 & YTemp!=BTemp])*log((BTemp[YTemp!=0 & YTemp!=BTemp]-YTemp[YTemp!=0 & YTemp!=BTemp])/(BTemp[YTemp!=0 & YTemp!=BTemp]-muTemp[YTemp!=0 & YTemp!=BTemp])))
                        DevianceRespTemp[YTemp==0]<-2*(BTemp[YTemp==0]*log((BTemp[YTemp==0])/(BTemp[YTemp==0]-muTemp[YTemp==0])))
                        DevianceRespTemp[YTemp==BTemp]<-2*(YTemp[YTemp==BTemp]*log(YTemp[YTemp==BTemp]/muTemp[YTemp==BTemp]))
                    }
                    if (RespDist[i]=="Poisson"){
                        DevianceRespTemp[YTemp!=0]<-2*(YTemp[YTemp!=0]*log(YTemp[YTemp!=0]/muTemp[YTemp!=0])-(YTemp[YTemp!=0]-muTemp[YTemp!=0]))
                        DevianceRespTemp[YTemp==0]<-2*muTemp[YTemp==0]
                    }
                    if (RespDist[i]=="Normal"){
                        DevianceRespTemp<-(YTemp-muTemp)^2
                    }
                    if (RespDist[i]=="Gamma"){
                        DevianceRespTemp<-2*(-log(YTemp/muTemp)+(YTemp-muTemp)/muTemp)
                        DiagPMATODEstim[(cModelsDims[i]+1):cModelsDims[i+1]]<-DiagPMATODEstim[(cModelsDims[i]+1):cModelsDims[i+1]]+1+as.vector(2*log(PhiTemp)/PhiTemp)+as.vector(2*digamma(1/PhiTemp)/PhiTemp)
                    }
                    DevianceResp[(cModelsDims[i]+1):cModelsDims[i+1]]<-DevianceRespTemp
                }
            }
            qrrO<-DiagPMAT[indexODEstim]
##            if (DEBUG) {print("qrrO");print(qrrO)}
            
            DevianceResp<-DevianceResp/(1-qrrO) 
##            if (DEBUG) {print("DevResp");print(DevianceResp)}
            # Algorithm for Gamma model #
            invSigmaGammaO<-((1-qrrO)/4)
        
            muGammaO<-Phi[indexODEstim]
        
            oldDYgamma<-DYgamma                
            # select DDY which are to estimate #
            DDYODEstim<-DDY[indexODEstim,,drop=F]
            # select columns which are not zero #
            tempind<-apply(matrix(as.logical(DDYODEstim),nrow(DDYODEstim),ncol(DDYODEstim)),2,any)
            indexgammaODEstim<-which(tempind==TRUE)
            # remove collums with all zeros from the design matrix #
            DDYODEstim<-DDYODEstim[,indexgammaODEstim,drop=F]
            DYgammaODEstim<-DYgamma[indexgammaODEstim]
            ksiO<-DDYODEstim%*%DYgammaODEstim
             
            ZOresp<-ksiO+(DevianceResp-muGammaO)/muGammaO
##            if (DEBUG) {print("ZOresp");print(ZOresp)}
            invSigmaGammaO<-diag(invSigmaGammaO)

            DYgammaODEstim<-solve(crossprod(DDYODEstim,invSigmaGammaO)%*%DDYODEstim,crossprod(DDYODEstim,invSigmaGammaO)%*%ZOresp)
            DYgamma[indexgammaODEstim]<-DYgammaODEstim
##            if (DEBUG) {print("DYgamma");print(DYgamma)}            
            Convergence<-Convergence+sum(abs(DYgamma-oldDYgamma))
##            if (DEBUG) {print("Convergence overdisp");print(Convergence)}
        }        
        
        VTCorr<-list(0)
        if (nrandcor>0){
            for (i in 1:length(CorrMat)){
                currentindex<-seq((cumqcorr[i]+1),cumqcorr[i+1])
                DRgamma[currentindex]<-log(diag(CholeskyMatrices[[i]]%*%t(CholeskyMatrices[[i]])))
                DRCorrgammaOld<-DRCorrgamma
                DRCorrgamma<-DRgamma[currentindex]   #HERE is probably an error if more than one correlation matrix is provided #
##                print("DRCorrgamma");print(DRCorrgamma)
                Convergence<-Convergence+sum(abs(DRCorrgammaOld-DRCorrgamma))
                # Transform the random effects back to the scale of correlated #
                VTemp<-unlist(VT)[(qcum[cumqcorr[i]+1]+1):qcum[cumqcorr[i+1]+1]] # Extract empirical bayes corresponding to the correlated effects of CorrMat[[i]]
                noraneff<-cumqcorr[i+1]-cumqcorr[i]
                VTemp<-matrix(VTemp,length(VTemp)/noraneff,noraneff)
                VTCorr[[i]]<-CholeskyMatrices[[i]]%*%t(VTemp)
                VTCorr[[i]]<-matrix(t(VTCorr[[i]]),nrow(VTemp)*noraneff,1)
                if (i==1) VTCorrTot<-VTCorr[[1]]
                else VTCorrTot<-c(VTCorrTot,VTCorr[[i]])
            }
            if (nrandind>0) VTCorrTot<-c(VTCorrTot,unlist(VT)[(qcum[cumqcorr[length(CorrMat)+1]]+1):qcum[length(qcum)]]) # This adds up independent empirical bayes to VTCorrTot
        } 
        
        # In the estimation of the correlation the original design matrix is used #
        ################################################################
        ##### Estimate Correlations  - Final Step in the Algorithm #####
        ################################################################
     
        # How will we approach this ? #
        # NOTES NOTES NOTES #
        # 1. We need to transform the design to the correlated desging #
        # 2. Compute first and second order derivatives #
        # 3. Update the correlations #
        # 4. Transform back to the independent design #
        # It is not that easy all this #
     
        # Iterate over correlations - for each compute the derivative #
        # double iterate over correlations for each pair compute hessian #
        
        # Strategy 1 - Use the previous variances values not the updated ones #
        # We update the whole matrix of correlations after we have variances and correlations updated #
     
        # Strategy 2 - use already here the updated matrix of variances - WHICH IS BETTER ? #
     
        # We try strategy 1 it seems faster to program but is it faster to work? #
     
       if (!is.null(CorrMat)) {
            ncorr<-length(unlist(Correlations))
            # How many random effects correlated there are - nrandcor  #
            ScoreCorr<-rep(0,ncorr)
            HessCorr<-matrix(0,ncorr,ncorr)
            dvhatdcorr<-list(0)
            dSigmadcorr<-list(0)        # Derivative of variance covariance matrix with respect to rho #
            # Note the derivative of Sigma with respect to two correlations (second order derivative) is equal to zero #
            # We need to remember that the change is made not for the correlation but on Fishers z transform -Inf/Inf #
            ZF <-list(0) #This is Fishers Z
            dCorrdZF<-list(0)
            OldCorrelations<-Correlations
            # Compute score equations #
            # The order in the program is as in the Correlations #
            # Iterated over the number of correlations not random effects which are correlated #
            
            # REMARK : It is wrong here in the way we treat correlations - it is not a vector but a list !#
            for (i in 1:length(CorrMat)) {
                # Compute current Fisher Z #
                # if (Correlations[i]==1) ZF[i]<--Inf
                # else if (Correlations[i]==-1) ZF[i]<-Inf               
                for (j in 1:length(Correlations[[i]])){
                    # Make ZF as a vector not as a list here #
                    ZF[[i]][j]<-0.5*log((1-Correlations[[i]][j])/(1+Correlations[[i]][j]))
                    dCorrdZF[[i]][j] <- - 4 * exp(2*ZF[[i]][j]) /((1+exp(2*ZF[[i]][i]))^2)
                }
            }
            # Compute derivative of Total Sigma matrix with respect to rho #
            # Here the creation of CorrMatDer #          
            
            for (i in 1:ncorr){ # This says which derivative     
                CurrentMatrix<-sum(cumqcorr<i)
                for (j in 1:length(CorrMat)){
                    currentindex<-i-sum(cumqcorr[cumqcorr<i])
                    if (j==CurrentMatrix) {
                        TempCorrDerMat<-SigmaMat[[j]]
                        if (Correlations[[j]][currentindex]!=0) TempCorrDerMat[CorrMat[[j]]==currentindex]<-SigmaMat[[j]][CorrMat[[j]]==currentindex]/(Correlations[[j]][currentindex])
                        else {
                            standardDeviations<-as.vector(sqrt(diag(SigmaMat[[j]])))
                            standardDeviations<-standardDeviations%*%t(standardDeviations)
                            TempCorrDerMat[CorrMat[[j]]==currentindex]<-standardDeviations[CorrMat[[j]]==currentindex]
                        }
                        TempCorrDerMat[CorrMat[[j]]!=currentindex]<-0
                        TempCorrDerMat<-TempCorrDerMat%x%diag(lcorr[j])
                        if (j==1) dSigmadcorr[[i]]<-TempCorrDerMat
                        else dSigmadcorr[[i]]<-dbind(dSigmadcorr[[i]],TempCorrDerMat[[j]])
                    }
                    else {
                        TempCorrDerMat<-SigmaMat[[j]]*0
                        TempCorrDerMat<-TempCorrDerMat%x%diag(lcorr[j])
                        if (j==1) dSigmadcorr[[i]]<-TempCorrDerMat
                        else dSigmadcorr[[i]]<-dbind(dSigmadcorr[[i]],TempCorrDerMat[[j]]) 
                    }                             

                }
            }
            for (i in 1:length(RandDist)){
                if (i==1) { dWRdUTot<-DWRDU(U[[i]],RandDist[i])
                            d2WRdU2Tot<-D2WRDU2(U[[i]],RandDist[i])
                            WRTot<-WR[[1]]
                            UTot<-U[[1]]
                }
                else {  dWRdUTot<-c(dWRdUTot,DWRDU(U[[i]],RandDist[i]))
                        d2WRdU2Tot<-c(d2WRdU2Tot,D2WRDU2(U[[i]],RandDist[i])) 
                        WRTot<-c(WRTot,WR[[i]])
                        UTot<-c(UTot,U[[i]])
                }    
            }
            
            # We need to add the zeros for uncorrelated random effects #
            if (nrandind>0) dSigmadcorr[[i]]<-dbind(dSigmadcorr[[i]],diag(qcum[length(qcum)]-qcum[cumqcorr[length(cumqcorr)]])*0)
 
            invSigmaTotComp<-invSigmaTot


            if ((qcum[cumqcorr[length(CorrMat)+1]+1]+1)-qcum[length(qcum)]<0) invSigmaTotComp<-dbind(invSigmaTot,diag(ISIGMAMvec[(qcum[cumqcorr[length(CorrMat)+1]]+1):qcum[length(qcum)]]))
            invTT2temp<-solve(t(ZOriginal)%*%(ZOriginal*as.vector(Wvec/Phi))+invSigmaTotComp)
            # within this iteration nest another iteration where hessian and second derivative of vhat is computed # 
            # Is this derivative for exponential family or only for normal distribution #
            for (i in 1:ncorr){ 
                dvhatdcorr[[i]]<-as.vector(invTT2temp%*%(invSigmaTotComp%*%(dSigmadcorr[[i]]%*%(invSigmaTotComp%*%VTCorrTot)))) # This seems to be okay
            }
            
            VTCorrTot<-as.vector(VTCorrTot)
            # DD matrix #
            MIDMAT<-dbind(diag(as.vector(Wvec/Phi)),invSigmaTot)
            if ((qcum[cumqcorr[length(CorrMat)+1]+1]+1)-qcum[length(qcum)]<0) MIDMAT<-dbind(MIDMAT,diag(ISIGMAMvec[(qcum[cumqcorr[length(CorrMat)+1]+1]+1):qcum[length(qcum)]]))
            DD<-t(TTOriginal)%*%MIDMAT%*%TTOriginal
 
            invDD<-solve(DD)
            
            # Computing score and hessian with respect to the correlation parameters #
            for (i in 1:ncorr) {
                #dDDdrho1<-dbind(matrix(0,ptot,ptot),-invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp)  ##### ERROR / FIXED - This implementation is only for normal random effects #####
                
                # The below version is for GLM distribution #
                dDDdrho1 <- dDDdranmat(X,ZOriginal,dWdmu,Wvec,dvhatdlambda=dvhatdcorr[[i]],invSigmaMat=invSigmaTotComp,dSigmadlambda=dSigmadcorr[[i]],WR=WRTot,dWRdu=dWRdUTot)
                
                # The below score is okay but must be adjusted for indepnedent random effects #
                ScoreCorr[i]<-(dvhatdcorr[[i]])%*%t(ZOriginal)%*%((Y-mu)*(1/Phi))-0.5*sum(diag(invSigmaTotComp%*%dSigmadcorr[[i]]))-dvhatdcorr[[i]]%*%invSigmaTotComp%*%VTCorrTot+
                    0.5*VTCorrTot%*%invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp%*%VTCorrTot-0.5*sum(diag(invDD%*%dDDdrho1)) # Score is okay ! #
              
                #ScoreCorr[i]<--0.5*sum(diag(invSigmaTotComp%*%dSigmadcorr[[i]]))+0.5*VTCorrTot%*%invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp%*%VTCorrTot-0.5*sum(diag(invDD%*%dDDdrho1))
              
                # Now proceed to the Hessian computation #
                for (j in 1:ncorr) {
                    #dDDdrho2<-dbind(matrix(0,ptot,ptot),-invSigmaTotComp%*%dSigmadcorr[[j]]%*%invSigmaTotComp)
                    dDDdrho2 <- dDDdranmat(X,ZOriginal,dWdmu,Wvec,dvhatdlambda=dvhatdcorr[[j]],invSigmaMat=invSigmaTotComp,dSigmadlambda=dSigmadcorr[[j]],WR=WRTot,dWRdu=dWRdUTot)
                    d2vhatdcorr2<--invTT2temp%*%((t(ZOriginal*as.vector(1/Phi)*dWdmu*Wvec*as.vector(ZOriginal%*%dvhatdcorr[[i]]))%*%as.vector(ZOriginal%*%dvhatdcorr[[j]]))+invSigmaTotComp%*%(dSigmadcorr[[j]]%*%(invSigmaTotComp%*%(dSigmadcorr[[i]]%*%(invSigmaTotComp%*%VTCorrTot))))+invSigmaTotComp%*%(dSigmadcorr[[i]]%*%
                            (invSigmaTotComp%*%(dSigmadcorr[[j]]%*%(invSigmaTotComp%*%VTCorrTot))))-invSigmaTotComp%*%(dSigmadcorr[[j]]%*%(invSigmaTotComp%*%dvhatdcorr[[i]]))-
                            invSigmaTotComp%*%(dSigmadcorr[[i]]%*%(invSigmaTotComp%*%dvhatdcorr[[j]])))
                    d2vhatdcorr2<-as.vector(d2vhatdcorr2)
                    #d2DDdrho2<-dbind(matrix(0,ptot,ptot),invSigmaTotComp%*%dSigmadcorr[[j]]%*%invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp+invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp%*%dSigmadcorr[[j]]%*%invSigmaTotComp)
                    
                    d2DDdrho2 <- d2DDdranmat2(X,ZOriginal,d2Wdmu2,dWdmu,Wvec,dvhatdlambda1=dvhatdcorr[[i]],dvhatdlambda2=dvhatdcorr[[j]],d2vhatdlambda12=d2vhatdcorr2,
                                invSigmaMat=invSigmaTotComp,dSigmadlambda1=dSigmadcorr[[i]],dSigmadlambda2=dSigmadcorr[[j]],d2Sigmadlambda12=dSigmadcorr[[j]]*0,WR=WRTot,dWRdu=dWRdUTot,d2WRdu2=d2WRdU2Tot)
                    HessCorr[i,j]<--dvhatdcorr[[i]]%*%t(ZOriginal)%*%(ZOriginal*as.vector(Wvec/Phi))%*%dvhatdcorr[[j]]+d2vhatdcorr2%*%t(ZOriginal*as.vector(1/Phi))%*%(Y-mu)+
                            0.5*sum(diag(invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp%*%dSigmadcorr[[j]]))+
                            dvhatdcorr[[i]]%*%invSigmaTotComp%*%dSigmadcorr[[j]]%*%invSigmaTotComp%*%VTCorrTot-dvhatdcorr[[j]]%*%invSigmaTotComp%*%dvhatdcorr[[i]]-d2vhatdcorr2%*%invSigmaTotComp%*%VTCorrTot+
                            dvhatdcorr[[j]]%*%invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp%*%VTCorrTot-
                            VTCorrTot%*%invSigmaTotComp%*%dSigmadcorr[[j]]%*%invSigmaTotComp%*%dSigmadcorr[[i]]%*%invSigmaTotComp%*%VTCorrTot-
                            # Adding the trace of the derivative #
                            0.5*(sum(diag(invDD%*%d2DDdrho2))-sum(diag(invDD%*%dDDdrho1%*%invDD%*%dDDdrho2)))
        
                
                }
        
            }
            # Now Lets create the derivative with respect to Fishers Z #
            dCorrdZF<-unlist(dCorrdZF)
            ScoreCorrZF<-ScoreCorr*dCorrdZF
            HessCorrZF<-matrix(0,ncorr,ncorr)
            for (i in 1:ncorr) {
                for (j in 1:ncorr) {
                    HessCorrZF[i,j]<-HessCorr[i,j]*dCorrdZF[i]*dCorrdZF[j]
                }
            }
            HessCorrZF<-(HessCorrZF+t(HessCorrZF))/2
            # Now we apply Newton-Raphson updater #
                        object1<-AdjProfCorrelations(unlist(ZF))
                        REML<-2*object1[[1]]
##                        print("REML");print(REML)
            if (EstimateCorrelations==TRUE){
                ZFOld<-ZF     
                if (nrow(as.matrix(-HessCorrZF))==1 & ncol(as.matrix(-HessCorrZF))==1) {
                    if (-HessCorrZF<=0) HessCorrZF<--0.00001    
                }
                else {
##                    HessCorrZF<--HessCorrZF
##                    print("HessCorrZF");print(HessCorrZF)
##                    print("ScoreCorrZF");print(ScoreCorrZF)
##                    if (det(-as.matrix(HessCorrZF))<=0) HessCorrZF<--nearPD(-as.matrix(HessCorrZF))       
                }
                check<-0
                frac<-1
                IHessCorrZF<-solve(HessCorrZF)
                object2<-AdjProfCorrelations(unlist(ZFOld))
                while (check==0) {    
                    TempZFStart<-unlist(ZF)-frac*(IHessCorrZF%*%ScoreCorrZF)
                    # If condition to check if the likelihood is increasing #
                
                    if (any(abs(TempZFStart)>3)) {
                        check<-0
                        frac<-frac/2
                    }
                    else {
                        object1<-AdjProfCorrelations(TempZFStart)
                        REML<-2*object1[[1]]
##                        print("REML");print(REML)
                        if (object1[[1]]-object2[[1]]<=0) {
                            check<-1
                            for (i in 1:length(CorrMat)){
                                ZF[[i]]<-TempZFStart[(1+cumqcorrels[i]):cumqcorrels[i+1]]   
                            }
                        }    
                        else {
                            check<-0
                            frac<-frac/2
                        }
                    }
                }
            ##### Update Correlations and design matrices #####
                M2h<--2*object1[[2]]
                M2pvh<--2*object1[[3]]
                M2pbvh<-REML
                CAIC<-object1[[4]]
                for (i in 1:length(CorrMat)){
                    Correlations[[i]]<-((1-exp(2*ZF[[i]]))/(1+exp(2*ZF[[i]])))
                }
            }
##            if (DEBUG==TRUE) print("Correlations");print(Correlations);
            Convergence<-Convergence+sum(abs(unlist(Correlations)-unlist(OldCorrelations)))
                M2h<--2*object1[[2]]
                M2pvh<--2*object1[[3]]
                M2pbvh<-REML
                CAIC<-object1[[4]]
            # Now make a re-design of the Z matrix taking into account new correlations and new variances #
            # We have new correlations - update correlations #
            
            TempCorrMat<-list(0)
            for (i in 1:length(CorrMat)){
                TempCorrMat[[i]]<-CorrMat[[i]]
                for (j in 1:length(Correlations[[i]])){
                    TempCorrMat[[i]][CorrMat[[i]]==j]<-Correlations[[i]][j]
                }
                diag(TempCorrMat[[i]])<-1
                CorrMatOut[[i]]<-TempCorrMat[[i]]
            }
   
            LambdaCorr<-exp(DDRCorr%*%DRCorrgamma)
            # now create correlation matrix for all the random effects #

            # create index of individual random effects #
       
            SigmaMat<-list(0)
            for (i in 1:length(CorrMat)){   
                LambdaLocal<-rep(0,qcorr[i])
                for (j in 1:qcorr[i]){
                    LambdaLocal[j]<-LambdaCorr[cumindCorrIndex[cumqcorr[i]+j]+1]
                }
                SigmaMat[[i]]<-sqrt(LambdaLocal)*t(CorrMatOut[[i]]*sqrt(LambdaLocal))
            }
            # merging the matrices by diagonals #
            for (i in 1:length(CorrMat)){
                if (i==1) {
                    SigmaTot<-SigmaMat[[1]]%x%diag(lcorr[1])
                    invSigmaTot<-solve(SigmaMat[[1]])%x%diag(lcorr[1])   
                }
                else {
                    SigmaTot<-dbind(SigmaTot,SigmaMat[[i]]%x%diag(lcorr[i]))
                    invSigmaTot<-dbind(invSigmaTot,solve(SigmaMat[[i]])%x%diag(lcorr[i]))
                }
            }
    
            # Matrix SigmaTot is the resulting matrix #
            # We have to make the random effects indendent via Cholesky decomposition #
            # We need a kroneker product of cholesky matrix times the SigmaTot #
            # You can do the cholesky on the total sigma matrix - the problem is the dimension is greater so maybe we loose computational efficiency #
            # DO cholesky on SigmaMat list #
            # We have to make functions to convert the ZZCorr into vectoral design for correlated random effects and back to the diagnoal design according to subject#
            # This will make the computationally more efficient things #
    
    
            # Now we modify the design matrix via cholesky decompositions #
            # All these steps need to be reprogramed using matrix package although in the final product there might be not so many zeros sometime #
            
            
            # REMARK: In what follows maybe not all steps are required to do again - if it takes a long time we can remove some - now i keep it for first attempt #
            for (i in 1:length(CorrMat)){
                itchol<-t(chol(SigmaMat[[i]]))  # This is actually cholesky decomposition instead of inverse, there was before inverse which was wrong
                CholeskyMatrices[[i]]<-itchol
                currentindex<-seq(cumqcorr[i]+1,cumqcorr[i+1])
                lengthcorrModelCur<-length(RespDist)
                valuescorrModelCur<-as.numeric(names(table(corrModel)))
                ZZCorrTemp<-rep(list(0),lengthcorrModelCur)
                for (j in 1:lengthcorrModelCur){
                    for (k in currentindex) {
                        if (ZZCorrTemp[[j]][1]==0 & length(ZZCorrTemp[[j]])==1){
                            if (corrModel[k]==valuescorrModelCur[j]) ZZCorrTemp[[j]]<-ZZCorrVec[[k]]
                            else ZZCorrTemp[[j]]<-rep(0,ModelsDims[j])
                        }
                        else {
                            if (corrModel[k]==valuescorrModelCur[j]) ZZCorrTemp[[j]]<-cbind(ZZCorrTemp[[j]],ZZCorrVec[[k]])
                            else ZZCorrTemp[[j]]<-cbind(ZZCorrTemp[[j]],rep(0,ModelsDims[j]))           
                        }
                    }
                }   
                # Binding it all together #
                for (j in 1:length(ZZCorrTemp)){
                    if (j==1) {
                        ZZCorrTempTot<-ZZCorrTemp[[j]]
                        nrowtot<-nrow(ZZCorrTemp[[j]])           
                    }
                    else {
                        ZZCorrTempTot<-rbind(ZZCorrTempTot,ZZCorrTemp[[j]])
                        nrowtot<-c(nrowtot,nrow(ZZCorrTemp[[j]]))
                    }
                }
                cnrowtot<-cumsum(c(0,nrowtot))
                # Now we use cholesky transform on the design matrix #
            
                ZZCorrTempTotUpd<-ZZCorrTempTot%*%itchol
                ZZShort[[i]]<-ZZCorrTempTot
                # ZZCorrTempTotUpd is the new design matrix for the joint model from the correlated part  #
                # This design matrix is in the short form (vector form) - we need to expand it to the diagnoal form # 
                # Now we need to take into accout to which model we should link the SSC #
        
        
                # Expansion to diagonal here we need matrix package already #
            
                for (j in currentindex){
                    DiagDesign[[j]]<-model.matrix(~SSC[[j]]-1)
                    DiagDesign[[j]]<-DiagDesign[[j]]*ZZCorrTempTotUpd[,(1+j-currentindex[1])]
                    if (j==currentindex[1]) ZZCorrTotUpd<-DiagDesign[[j]]
                    else ZZCorrTotUpd<-cbind(ZZCorrTotUpd,DiagDesign[[j]])
                    #DiagDesign[[j]]<-model.matrix(~SSC[[j]]-1)
                }
        
                ZZCorrUpd[[i]]<-ZZCorrTotUpd 
        
                if (i==1) ZZCorrDesign<-ZZCorrUpd[[1]]
                else ZZCorrDesign<-cbind(ZZCorrDesign,ZZCorrUpd[[i]])
                # Fill with zeros the ZZCorrDesign #
                ZZCorrDesign<-FillZeros(ZZCorrDesign,nr=sum(ModelsDims),nc=ncol(ZZCorrDesign))
            }
            q<-sapply(DiagDesign,ncol)
        # End of if corrmat not null #
        }
        
        if (nrandind>0) {
            ZZIndepTemp<-list(0)
            for (i in 1:nrandind){
                if (indepModel[i]==length(ModelsDims)) rowsdown<-0
                else rowsdown<-sum(ModelsDims[(indepModel[i]+1):length(ModelsDims)])
                if (indepModel[i]==1) rowsup<-0
                else rowsup<-sum(ModelsDims[1:(indepModel[i]-1)])
            
                # First expand up if possible #
                ZZIndepTemp[[i]]<-FillZeros(ZZIndep[[i]],nr=(nrow(ZZIndep[[i]])+rowsup),nc=ncol(ZZIndep[[i]]),up=TRUE)
                ZZIndepTemp[[i]]<-FillZeros(ZZIndepTemp[[i]],nr=(nrow(ZZIndepTemp[[i]])+rowsdown),nc=ncol(ZZIndepTemp[[i]]))
                if (!exists("ZZCorrDesign") & i==1) ZZDesign<-ZZIndepTemp[[i]]
                if (exists("ZZCorrDesign") & i==1) ZZDesign<-cbind(ZZCorrDesign,ZZIndepTemp[[1]])
                if (i>1) ZZDesign<-cbind(ZZDesign,ZZIndepTemp[[i]])
            }
        
        }
        else if(nrandcor>0) ZZDesign<-ZZCorrDesign
        TT<-cbind(X,ZZDesign)
        if (nrandcor>0){
            for (i in 1:nrandcor){
                TT<-rbind(TT,cbind(matrix(0,q[i],ptot+qcum[i]),diag(q[i]),matrix(0,q[i],qcum[nrandcor+nrandind+1]-qcum[i+1])))
            }
        }
    
        if (nrandind>0){
            for (i in 1:nrandind){
                TT<-rbind(TT,cbind(matrix(0,q[i+nrandcor],ptot+qcum[i+nrandcor]),diag(q[i+nrandcor]),matrix(0,q[i+nrandcor],qcum[nrandind+nrandcor+1]-qcum[nrandcor+i+1])))
            }
        }  
    
        Z<-TT[(1:sum(ModelsDims)),((ptot+1):ncol(TT))]                      
            
    }        
    if (StandardErrors==TRUE) {
    # Now we need standard errors and diagnostics and we are done #
    for (i in 1:length(RandDist)){
        if (i==1) { 
            if (nrandcor>0) U[[i]]<-LinkR(VTCorrTot[(qcum[i]+1):qcum[i+1]],RandDist[i])
            dWRdUTot<-DWRDU(U[[i]],RandDist[i])
            d2WRdU2Tot<-D2WRDU2(U[[i]],RandDist[i])
            WRTot<-WR[[1]]
            UTot<-U[[1]]
            bfuncv<- bfuncvgen(Vvec=V[[1]],Dist=RandDist[i])
        }
        else {  
            if (nrandcor>0) U[[i]]<-LinkR(VTCorrTot[(qcum[i]+1):qcum[i+1]],RandDist[i])
            dWRdUTot<-c(dWRdUTot,DWRDU(U[[i]],RandDist[i]))
            d2WRdU2Tot<-c(d2WRdU2Tot,D2WRDU2(U[[i]],RandDist[i])) 
            WRTot<-c(WRTot,WR[[i]])
            UTot<-c(UTot,U[[i]])
            bfuncv<-c(bfuncv,bfuncvgen(Vvec=V[[i]],Dist=RandDist[i]))
        }    
    }
    # Standard errors gradient hessian - correlation parameters #
    # But these should be jointly with DRgamma #
    if (!is.null(CorrMat)){
        StdErrCorr<-sqrt(diag(solve(-HessCorr)))            # These dont take into account the variablity of the dispersion components #
        StdErrCorrZF<-sqrt(diag(solve(-HessCorrZF)))        # Correlation need to be jointly with DRgamma values and then inverted
    }
    else {
        StdErrCorr<-NULL
        StdErrCorrZF<-NULL
        invSigmaTotComp<-diag(ISIGMAMvec[(ntot+1):(ntot+qcum[nrandcor+nrandind+1])])
        invTT2temp<-solve(t(ZOriginal*as.vector(Wvec/Phi))%*%(ZOriginal)+invSigmaTotComp)
    }
        
    # Standard errors for random effects - these are also not all ok should be jointly with fixed effects #
    
    
    
    
    # Standard errors for Beta #
    # Now we have to consider two situations: one when LAPFIX is used and one when it is not #
    # It is now a problem if some Betas are estimated without the Laplace Approximations and some are not #
    # clearly the Beta estimates are not orthogonal among each other, therefore for now we should adapt the approach that either 
    # all betas are estimated with the adjusted profile likelihood or with a joint likelihood #
    
    if (all(!LaplaceFixed)==TRUE) {
        StdErrBeta<-sqrt(diag(invDD)[1:ptot])
        StdErrVs<-sqrt(diag(invDD)[(ptot+1):length(diag(invDD))]) # These standard errors are not the same as standard errors of PROC MIXED, the part of them is the same!
        
    }                    
    
    else if (all(LaplaceFixed)==TRUE) {
        # Compute the derivative of vhat with respect to beta - it is a 3 dimensional array q (number of random effects) times (p x p) #
        # To musi byc zmienione jesli jest wiecej niz jeden model !!!!!! #
         
        dvhatdbeta<--invTT2temp%*%(t(ZOriginal*as.vector(Wvec/Phi)))%*%X   # Pytanie czy ta pochodna poprawnie traktuje vhat z modelu (1) wzgledem beta z modelu (2) czy sa to zera? ale czy powinny to byc zera? 
                                                                                # Pochodna jest OK prawdopodobnie
        # !!!!!! We need to use invSigmaTotComp when we use normal random effects and correlation - but is it multiplied by WR the independent parts of random effects ? !!!!! #
        # Yes invSigmaTotComp is dbind with ISIGMAMvec, which is the multiplication #
        # Compute second derivatives of dvhatdbeta1dbeta2 #
        # You have to remember that now the matrix is from correlation of random effects #
        # @!!!!! DWRDU is not stored anywhere you need to use a function for that #
        # Also WR check where this is stored #
        d2vhatdbeta2<-array(0,dim=c(ptot,ptot,length(VT)))
        # compute dWRdUTot #
            for (i in 1:ptot) {
                for (j in 1:ptot) {
                # This is more complicated when more than one model is used  #
                # Proper adjustments need to be made #
                # This diagonal is tough to compute #
                
                # IF both beta from the same model  then the derivative is ok #
                
                #DiagDesignMat <- 1 ###                      # Macierz ktora jest Xp+Zdvhatdbeta - definiuje ktore rzeczy sa na diagonals - 
                d2vhatdbeta2[i,j,]<--invTT2temp%*%(t(ZOriginal*as.vector((dWdmu*Wvec)/Phi))%*%((X[,i]+ZOriginal%*%dvhatdbeta[,i])*(X[,j]+ZOriginal%*%dvhatdbeta[,j]))+
                                                diag(as.vector((dWRdUTot*WRTot)/Lambda))%*%(dvhatdbeta[,i]*dvhatdbeta[,j])) # This line here needs to be ammended #
                #d2vhatdbeta2[i,j,]<--invTT2temp%*%t(ZOriginal*as.vector((dWdmu*Wvec)/Phi)%*% # Nowa definicja 
                }
            }    
           # end if nModels equal 1
        # Below we implement the computation of derivatives for vhat for more than one model !!!!
             
        
            # Now we compute the hessian of the Laplace Approximation to the marginal likelihood #
        
            # We need hessian of the TT2 matrix with respect to Beta and Beta*Beta #
        
            # For the hessian we need the second derivative of the Wvec matrix #
        
            # We directly compute the hessian of the Laplace Approximation #
        HessianLapFixBeta<-matrix(0,ptot,ptot)
        for (i in 1:ptot) {
            for (j in 1:i) {
                
                dTT2dbetai<-t(ZOriginal*as.vector((1/Phi)*dWdmu*Wvec*(X[,i]+ZOriginal%*%dvhatdbeta[,i])))%*%ZOriginal+diag(as.vector((1/Lambda)*dWRdUTot*WRTot*dvhatdbeta[,i]))
                dTT2dbetaj<-t(ZOriginal*as.vector((1/Phi)*dWdmu*Wvec*(X[,j]+ZOriginal%*%dvhatdbeta[,j])))%*%ZOriginal+diag(as.vector((1/Lambda)*dWRdUTot*WRTot*dvhatdbeta[,j]))
                d2TT2dbetaidbetaj<-t(ZOriginal*as.vector((1/Phi)*d2Wdmu2*Wvec*(X[,i]+ZOriginal%*%dvhatdbeta[,i])*Wvec*(X[,j]+ZOriginal%*%dvhatdbeta[,j])))%*%ZOriginal+
                                        t(ZOriginal*as.vector((1/Phi)*dWdmu*dWdmu*Wvec*(X[,i]+ZOriginal%*%dvhatdbeta[,i])*(X[,j]+ZOriginal%*%dvhatdbeta[,j])))%*%ZOriginal+
                                        t(ZOriginal*as.vector((1/Phi)*dWdmu*Wvec*(ZOriginal%*%d2vhatdbeta2[i,j,])))%*%ZOriginal+
                                        diag(as.vector((1/Lambda)*d2WRdU2Tot*WRTot*dvhatdbeta[,i]*WRTot*dvhatdbeta[,j]))+
                                        diag(as.vector((1/Lambda)*dWRdUTot*dWRdUTot*WRTot*dvhatdbeta[,i]*dvhatdbeta[,j]))+
                                        diag(as.vector((1/Lambda)*dWRdUTot*WRTot*d2vhatdbeta2[i,j,]))
                # !!!!!!! In the above we need to replace Lambda by invSigmaTotComp - although it does not matter as for normal dWRduTot=0 and dWdmu=0 #
                
                # Maybe in Hessian Lap Fix we can also add terms which sum up to zero for better accuracy ?#
                HessianLapFixBeta[i,j]<-sum(-(as.vector((X[,i]+ZOriginal%*%dvhatdbeta[,i])*Wvec*(1/Phi)*(X[,j]+ZOriginal%*%dvhatdbeta[,j])))) + sum(as.vector((1/Phi)*(Y-mu)*(ZOriginal%*%d2vhatdbeta2[i,j,])))-
                                                as.vector(dvhatdbeta[,i])%*%invSigmaTotComp%*%as.vector(dvhatdbeta[,j])+as.vector((PsiM-UTot))%*%(invSigmaTotComp)%*%as.vector(d2vhatdbeta2[i,j,])-
                                                0.5*sum(diag(invTT2temp%*%d2TT2dbetaidbetaj))+0.5*sum(diag(invTT2temp%*%dTT2dbetai%*%invTT2temp%*%dTT2dbetaj))
                                                
                                                
                HessianLapFixBeta[j,i]<-HessianLapFixBeta[i,j]
            }
        }
               
        StdErrBeta<-sqrt(diag(solve(-HessianLapFixBeta)))   
        StdErrVs<-sqrt(diag(invTT2temp))    
    }
    
    else {
        warning("Some fixed effects are estimated by Laplace some by h-likelihood; currently we do not know how to compute standard errors in this case")
        StdErrBeta<-NULL
    } 
    
    # Now we compute the standard errors of residual dispersion components #
    
    # The problem is to determine which components are estimated and which are kept fixed #
    SelectGamma<-NULL
    SelectModel<-NULL
    StdErrODEst <- NULL
    for (i in 1:nModels) {
        if (EstimateOverDisp[i] == TRUE) {
            rowsTemp <- (cModelsDims[i]+1):(cModelsDims[i+1])
            TempOD <- as.matrix(DDY[rowsTemp,])
            colnames(TempOD)<-1:ncol(TempOD)
            columnsTemp <- apply(matrix(as.logical(TempOD),nrow(TempOD),ncol(TempOD)),2,any)
            whichGamma <- colnames(TempOD)[which(columnsTemp==TRUE)]
            SelectModelTemp <- rep(i,length(whichGamma))
            if (i==1) {
                SelectGamma <- as.numeric(whichGamma)
                SelectModel <- SelectModelTemp
            }
            else {
                SelectGamma <- c(SelectGamma, as.numeric(whichGamma))
                SelectModel <- c(SelectModel, SelectModelTemp)   
            }     # This says which gamma are estimated and with respect to them the hessian is going to be computed #
        }  
    }
    HessianODEst<-matrix(0,length(SelectGamma),length(SelectGamma))
    
    DevianceRespTotal <- rep(0,ntot)  
    for (i in 1:nModels) {
        DevianceRespTemp<-rep(0,cModelsDims[i+1]-cModelsDims[i])
        YTemp<-YList[[i]]
        BTemp<-B[(cModelsDims[i]+1):cModelsDims[i+1]]
        muTemp<-mu[(cModelsDims[i]+1):cModelsDims[i+1]]
        PhiTemp<-Phi[(cModelsDims[i]+1):cModelsDims[i+1]]
        if (RespDist[i]=="Binomial") {
            DevianceRespTemp[YTemp!=0 & YTemp!=BTemp]<-2*(YTemp[YTemp!=0 & YTemp!=BTemp]*log(YTemp[YTemp!=0 & YTemp!=BTemp]/muTemp[YTemp!=0 & YTemp!=BTemp])-(YTemp[YTemp!=0 & YTemp!=BTemp]-BTemp[YTemp!=0 & YTemp!=BTemp])*log((BTemp[YTemp!=0 & YTemp!=BTemp]-YTemp[YTemp!=0 & YTemp!=BTemp])/(BTemp[YTemp!=0 & YTemp!=BTemp]-muTemp[YTemp!=0 & YTemp!=BTemp])))
            DevianceRespTemp[YTemp==0]<-2*(BTemp[YTemp==0]*log((BTemp[YTemp==0])/(BTemp[YTemp==0]-muTemp[YTemp==0])))
            DevianceRespTemp[YTemp==BTemp]<-2*(YTemp[YTemp==BTemp]*log(YTemp[YTemp==BTemp]/muTemp[YTemp==BTemp]))
        }
        if (RespDist[i]=="Poisson"){
            DevianceRespTemp[YTemp!=0]<-2*(YTemp[YTemp!=0]*log(YTemp[YTemp!=0]/muTemp[YTemp!=0])-(YTemp[YTemp!=0]-muTemp[YTemp!=0]))
            DevianceRespTemp[YTemp==0]<-2*muTemp[YTemp==0]
        }
        if (RespDist[i]=="Normal"){
            DevianceRespTemp<-(YTemp-muTemp)^2
        }
        if (RespDist[i]=="Gamma"){
            DevianceRespTemp<-2*(-log(YTemp/muTemp)+(YTemp-muTemp)/muTemp)
            #DiagPMATODEstim[(cModelsDims[i]+1):cModelsDims[i+1]]<-DiagPMATODEstim[(cModelsDims[i]+1):cModelsDims[i+1]]+1+as.vector(2*log(PhiTemp)/PhiTemp)+as.vector(2*digamma(1/PhiTemp)/PhiTemp)
        }
        if (i == 1) DevianceRespTotal <- DevianceRespTemp
        else DevianceRespTotal <- c(DevianceRespTotal, DevianceRespTemp)
    }
    
    # The part below needs to be removed #
    #DevD<-DevianceRespTotal
    #PhiD<-Phi
    #XD<-X
    #ZD<-ZOriginal
    #WD<-Wvec
    #invD<-invSigmaTotComp
    
    DDr1<-cbind(t(X*as.vector(Wvec/Phi))%*%X,t(X*as.vector(Wvec/Phi))%*%ZOriginal)
    DDr2<-cbind(t(ZOriginal*as.vector(Wvec/Phi))%*%X,t(ZOriginal*as.vector(Wvec/Phi))%*%ZOriginal+invSigmaTotComp)
    DD<-rbind(DDr1,DDr2)
    solveDD<-solve(DD)
    # Now compute the actual Hessian #
    if (!is.null(SelectGamma)) {

        for (i in 1:length(SelectGamma)) {
            for (j in 1:i){
                modelPhi1<-SelectModel[i]
                modelPhi2<-SelectModel[j]
                if (modelPhi1 != modelPhi2 ) { 
                    d2Qdphi2 <- 0
                    d2Qdgamma2 <- 0
                }
                # Derivative of the Quasi likelihood #
                else {
                    d2Qdphi2 <- -(DevianceRespTotal[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]]/(Phi[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]]^3))+
                                    0.5*(1/(Phi[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]]^2))   
                     
                    if (RespDist[SelectModel[i]]=="Gamma") {
                        PhiGCur<-as.vector(Phi[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]])
                        d2Qdphi2<-d2Qdphi2-(2*log(PhiGCur))/(PhiGCur^3)+(1/PhiGCur^3)-0.5*(1/PhiGCur^2)-(2*digamma(1/PhiGCur))/(PhiGCur^3)-(trigamma(1/PhiGCur)/PhiGCur^4)
                    }
                    #DDYODCuri<-rep(0,ntot)
                    DDYODCuri<-DDY[(cModelsDims[SelectModel[i]]+1):(cModelsDims[SelectModel[i]+1]),SelectGamma[i]]
                    #DDYODCurj<-rep(0,ntot)
                    DDYODCurj<-DDY[(cModelsDims[SelectModel[j]]+1):(cModelsDims[SelectModel[j]+1]),SelectGamma[j]]                
                
                    d2Qdgamma2 <- (d2Qdphi2 * (as.vector(Phi[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]])^2) * as.vector(DDYODCuri))%*%
                            as.vector(DDYODCurj)
##                    cat("i: ",i,"; j:",j," s2Qdgamma2:",d2Qdgamma2)
                }
                            
                # !!!!! Here a correction for the gamma distribution so the h is used instead of Q !!!!!! #
                            
                # Derivative of the determinant #
                if (modelPhi1 == modelPhi2) {
                    PhiCur <- rep (0,ntot+qtot)
                    PhiCur[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]]<-Phi[(cModelsDims[SelectModel[i]]+1):cModelsDims[SelectModel[i]+1]]
                    dimmOD<-nrow(TTOriginal)-length(Wvec)
                    WvecODCur<-c(Wvec,rep(0,dimmOD))
                    
                    DDYODCuri<-rep(0,ntot)
                    DDYODCuri[(cModelsDims[i]+1):(cModelsDims[i+1])]<-DDY[(cModelsDims[SelectModel[i]]+1):(cModelsDims[SelectModel[i]+1]),SelectGamma[i]]
                    DDYODCurj<-rep(0,ntot)
                    DDYODCurj[(cModelsDims[j]+1):(cModelsDims[j+1])]<-DDY[(cModelsDims[SelectModel[j]]+1):(cModelsDims[SelectModel[j]+1]),SelectGamma[j]]                      
                    
                    dDDdgamma1p1 <- - cbind(t(X*as.vector(Wvec*DDYODCuri/Phi))%*%X,t(X*as.vector(Wvec*DDYODCuri/Phi))%*%ZOriginal)
                    dDDdgamma1p2 <- - cbind(t(ZOriginal*as.vector(Wvec*DDYODCuri/Phi))%*%X,t(ZOriginal*as.vector(Wvec*DDYODCuri/Phi))%*%ZOriginal)
                    dDDdgamma1    <-  rbind(dDDdgamma1p1,dDDdgamma1p2)
                    
                    dDDdgamma2p1 <- - cbind(t(X*as.vector(Wvec*DDYODCurj/Phi))%*%X,t(X*as.vector(Wvec*DDYODCurj/Phi))%*%ZOriginal)
                    dDDdgamma2p2 <- - cbind(t(ZOriginal*as.vector(Wvec*DDYODCurj/Phi))%*%X,t(ZOriginal*as.vector(Wvec*DDYODCurj/Phi))%*%ZOriginal)
                    dDDdgamma2   <-   rbind(dDDdgamma2p1,dDDdgamma2p2)
                    
                    d2DDdgamma12p1 <-  2* cbind(t(X*as.vector(Wvec*DDYODCuri*DDYODCurj/Phi))%*%X,t(X*as.vector(Wvec*DDYODCuri*DDYODCurj/Phi))%*%ZOriginal)
                    d2DDdgamma12p2 <-  2* cbind(t(ZOriginal*as.vector(Wvec*DDYODCuri*DDYODCurj/Phi))%*%X,t(ZOriginal*as.vector(Wvec*DDYODCuri*DDYODCurj/Phi))%*%ZOriginal)                    
                    d2DDdgamma12   <-  rbind(d2DDdgamma12p1,d2DDdgamma12p2)
                    
                    HessianODEst[i,j]<-d2Qdgamma2-0.5*sum(diag(solveDD%*%d2DDdgamma12))+0.5*sum(diag(solveDD%*%dDDdgamma1%*%solveDD%*%dDDdgamma2))
                    HessianODEst[j,i]<-HessianODEst[i,j]
##                    cat("Parameter: ",i," Drugi: ",j,"\n d2Q: ",d2Qdgamma2," trace: ",-0.5*sum(diag(solveDD%*%d2DDdgamma12))+0.5*sum(diag(solveDD%*%dDDdgamma1%*%solveDD%*%dDDdgamma2)))
                }
                else {
                    DDYODCuri<-rep(0,ntot)
                    DDYODCuri[(cModelsDims[i]+1):(cModelsDims[i+1])]<-DDY[(cModelsDims[SelectModel[i]]+1):(cModelsDims[SelectModel[i]+1]),SelectGamma[i]]
                    DDYODCurj<-rep(0,ntot)
                    DDYODCurj[(cModelsDims[j]+1):(cModelsDims[j+1])]<-DDY[(cModelsDims[SelectModel[j]]+1):(cModelsDims[SelectModel[j]+1]),SelectGamma[j]]                      
                    
                    dDDdgamma1p1 <- - cbind(t(X*as.vector(Wvec*DDYODCuri/Phi))%*%X,t(X*as.vector(Wvec*DDYODCuri/Phi))%*%ZOriginal)
                    dDDdgamma1p2 <- - cbind(t(ZOriginal*as.vector(Wvec*DDYODCuri/Phi))%*%X,t(ZOriginal*as.vector(Wvec*DDYODCuri/Phi))%*%ZOriginal)
                    dDDdgamma1    <-  rbind(dDDdgamma1p1,dDDdgamma1p2)
                    
                    dDDdgamma2p1 <- - cbind(t(X*as.vector(Wvec*DDYODCurj/Phi))%*%X,t(X*as.vector(Wvec*DDYODCurj/Phi))%*%ZOriginal)
                    dDDdgamma2p2 <- - cbind(t(ZOriginal*as.vector(Wvec*DDYODCurj/Phi))%*%X,t(ZOriginal*as.vector(Wvec*DDYODCurj/Phi))%*%ZOriginal)
                    dDDdgamma2   <-   rbind(dDDdgamma2p1,dDDdgamma2p2)
                                   
                    HessianODEst[i,j]<-0.5*sum(diag(solveDD%*%dDDdgamma1%*%solveDD%*%dDDdgamma2))
                    HessianODEst[j,i]<-HessianODEst[i,j]
                }
            }    
        }
        colnames(HessianODEst)<-SelectGamma
        rownames(HessianODEst)<-SelectGamma
        StdErrODEst<-diag(sqrt(solve(-HessianODEst)))
        names(StdErrODEst)<-paste("gamma",SelectGamma,sep="")
    }
    
    ##########################################
    ##### Hessian for DRgamma parameters #####
    ##########################################
    

    # We can use DVhatDlambda as it is already calculated 
    if (nrandcor > 0) nrandtot<-nrandcor
    else nrandtot <- 0
    if (nrandind > 0) nrandtot<-nrandtot+nrandind
    
    HessianRVC<-matrix(0,nrandtot,nrandtot)
    
    # Hessian of correlated part #
    # debug(dhdranmat)
    # debug(dDDdranmat)
    # VTCorrTot as correlated random effects #
    
    if (nrandcor > 0) {
        DSigmadlambdaConstant<-solve(invSigmaTotComp)   # Copy current matrix #
        for (i in 1:length(CorrMat)){ 
            # Compute the derivative of dSigmadlambda #
            # Determine which lambda #
            LambdaLocal<-rep(0,qcorr[i])
            for (j in 1:qcorr[i]){    
                LambdaLocal[j]<-LambdaCorr[cumindCorrIndex[cumqcorr[i]+j]+1]
            }
            FDER <- rep(0,qcorr[i]) # gradient 
            for (Dindex in 1:qcorr[i]){ # says which derivative 
                CorrMatOutDeriv<-CorrMatOut[[i]]
                diag(CorrMatOutDeriv)<-rep(0,qcorr[i])
                CorrMatOutDeriv[Dindex,]<-CorrMatOutDeriv[Dindex,]/2
                CorrMatOutDeriv[,Dindex]<-CorrMatOutDeriv[,Dindex]/2
                CorrMatOutDeriv[Dindex,Dindex]<-1
                CorrMatOutDeriv[-Dindex,-Dindex]<-0
                DSigmaMatdlambda<-t(t(sqrt(LambdaLocal)*CorrMatOutDeriv)*(sqrt(LambdaLocal)/(LambdaLocal[Dindex])))
                DSigmaTotdlambda<-DSigmaMatdlambda%x%diag(lcorr[i])
                # Computing first order derivative #
                DSigmadlambda1<-matrix(0,qcum[nrandtot+1],qcum[nrandtot+1])
                DSigmadlambda1[(cumindCorrIndex[cumqcorr[i]+1]+1):(cumindCorrIndex[cumqcorr[i+1]+1]),(cumindCorrIndex[cumqcorr[i]+1]+1):(cumindCorrIndex[cumqcorr[i+1]+1])]<-
                    DSigmaTotdlambda
                dvhatdlambda1 <- as.vector(dvhatdranmat(invTT2=invTT2temp,invSigmaMat=invSigmaTotComp,dSigmadlambda=DSigmadlambda1,Psi=PsiM,Uvec=UTot))
                dDDdlambda1<-dDDdranmat(X=X,Z=ZOriginal,dWdmu=dWdmu,Wvec=Wvec,dvhatdlambda=dvhatdlambda1,invSigmaMat=invSigmaTotComp,dSigmadlambda=DSigmadlambda1,
                                    WR=WRTot,dWRdu=dWRdUTot)
                LambdaCur1<-LambdaLocal[Dindex]
                FDER[Dindex] <-dhdranmatCorr(Z=ZOriginal,y=Y,mu=mu,Phi=Phi,dvhatdlambda=dvhatdlambda1,invSigmaMat=invSigmaTotComp,Psi=PsiM,Uvec=UTot,
                                                Vvec=VTCorrTot,bfuncv=bfuncv,dSigmadlambda=DSigmadlambda1)-0.5*sum(diag(solveDD%*%dDDdlambda1)) 
                if (Dindex==1) {
                    dVcur1 <- dvhatdlambda1 
                    dSigm1 <- DSigmadlambda1
                    dDDd1 <- dDDdlambda1
                }
                # Now computing the second derivative #
                for (i2 in 1:length(CorrMat)) {
                    LambdaLocal<-rep(0,qcorr[i])
                    for (j in 1:qcorr[i2]){    
                        LambdaLocal[j]<-LambdaCorr[cumindCorrIndex[cumqcorr[i2]+j]+1]
                    }                
                    
                    for (Dindex2 in 1:qcorr[i2]){ # says which derivative 
                        CorrMatOutDeriv<-CorrMatOut[[i2]]
                        diag(CorrMatOutDeriv)<-rep(0,qcorr[i2])
                        CorrMatOutDeriv[Dindex2,]<-CorrMatOutDeriv[Dindex2,]/2
                        CorrMatOutDeriv[,Dindex2]<-CorrMatOutDeriv[,Dindex2]/2
                        CorrMatOutDeriv[Dindex2,Dindex2]<-1
                        CorrMatOutDeriv[-Dindex2,-Dindex2]<-0
                        DSigmaMatdlambda<-t(t(sqrt(LambdaLocal)*CorrMatOutDeriv)*(sqrt(LambdaLocal)/(LambdaLocal[Dindex2])))
                        DSigmaTotdlambda<-DSigmaMatdlambda%x%diag(lcorr[i2])
                        # Computing second first order derivative #
                        DSigmadlambda2<-matrix(0,qcum[nrandtot+1],qcum[nrandtot+1])
                        DSigmadlambda2[(cumindCorrIndex[cumqcorr[i2]+1]+1):(cumindCorrIndex[cumqcorr[i2+1]+1]),(cumindCorrIndex[cumqcorr[i2]+1]+1):(cumindCorrIndex[cumqcorr[i2+1]+1])]<-
                            DSigmaTotdlambda
                        dvhatdlambda2 <- as.vector(dvhatdranmat(invTT2=invTT2temp,invSigmaMat=invSigmaTotComp,dSigmadlambda=DSigmadlambda2,Psi=PsiM,Uvec=UTot))
                        
                        dDDdlambda2<-dDDdranmat(X=X,Z=ZOriginal,dWdmu=dWdmu,Wvec=Wvec,dvhatdlambda=dvhatdlambda2,invSigmaMat=invSigmaTotComp,dSigmadlambda=DSigmadlambda2,
                                    WR=WRTot,dWRdu=dWRdUTot)
                        LambdaCur2<-LambdaLocal[Dindex2]
                        # Computing the second order derivative of the SigmaMat over the lambda #
                        # Three possiblities: 1 / same matrix same sigma
                        #                   : 2 / same matrix different sigma
                        #                   : 3 / different matrix CorrMat
                        D2Sigmadlambda12<-matrix(0,qcum[nrandtot+1],qcum[nrandtot+1])
                        if (i==i2 & Dindex==Dindex2) {
                            D2SigmaMatdlambda12 <- DSigmaMatdlambda/(-2*LambdaLocal[Dindex2])
                            D2SigmaMatdlambda12[Dindex,Dindex]<-0
                            D2SigmaTotdlambda12 <- D2SigmaMatdlambda12%x%diag(lcorr[i2])
                            D2Sigmadlambda12[(cumindCorrIndex[cumqcorr[i2]+1]+1):(cumindCorrIndex[cumqcorr[i2+1]+1]),(cumindCorrIndex[cumqcorr[i2]+1]+1):(cumindCorrIndex[cumqcorr[i2+1]+1])]<-
                                D2SigmaTotdlambda12
                        }
                        if (i==i2 & Dindex!=Dindex2) {
                             D2SigmaMatdlambda12temp <- DSigmadlambda2/(2*LambdaLocal[Dindex])
                             D2SigmaMatdlambda12 <- matrix(0,qcorr[i],qcorr[i])
                             D2SigmaMatdlambda12[Dindex,Dindex2] <- D2SigmaMatdlambda12temp[Dindex,Dindex2]
                             D2SigmaMatdlambda12[Dindex2,Dindex] <- D2SigmaMatdlambda12temp[Dindex2,Dindex]   
                             D2SigmaTotdlambda12 <- D2SigmaMatdlambda12%x%diag(lcorr[i2])
                             D2Sigmadlambda12[(cumindCorrIndex[cumqcorr[i2]+1]+1):(cumindCorrIndex[cumqcorr[i2+1]+1]),(cumindCorrIndex[cumqcorr[i2]+1]+1):(cumindCorrIndex[cumqcorr[i2+1]+1])]<-
                                D2SigmaTotdlambda12
                              
                        }
                        if (i!=i2) {
                            D2SigmaMatdlambda12 <- matrix(0, qcorr[i], qcorr[i])
                            D2Sigmadlambda12 <- matrix(0, qcum[nrandtot+1], qcum[nrandtot+1])
                        }
                        
                        d2vhatdlambda12<-as.vector(d2vhatdranmat2(invTT2=invTT2temp,Z=ZOriginal,Phi=Phi,dWdmu=dWdmu,Wvec=Wvec,
                                dvhatdlambda1=dvhatdlambda1,dvhatdlambda2=dvhatdlambda2,invSigmaMat=invSigmaTotComp,WR=WRTot,
                                dWRdu=dWRdUTot,dSigmadlambda1=DSigmadlambda1,dSigmadlambda2=DSigmadlambda2,Psi=PsiM,Uvec=UTot,d2Sigmadlambda12=D2Sigmadlambda12))
                        
                        if (Dindex2==1) {
                            dVcur2<-dvhatdlambda2
                            dSigm2<-DSigmadlambda2
                            dDDd2<-dDDdlambda2
                        }

                        # Computing which derivative we deal with #
                        firstindex<-cumqcorr[i]+Dindex
                        secondindex<-cumqcorr[i2]+Dindex2
                        
                        # Computing the actual hessian #
                        
                        # Define 
                        d2hdlambda2 <- d2hdranmatCorrCorr(Z=ZOriginal,y=Y,mu=mu,Phi=Phi,d2vhatdlambda12=d2vhatdlambda12,dvhatdlambda1=dvhatdlambda1,
                            dvhatdlambda2=dvhatdlambda2,Wvec=Wvec,invSigmaMat=invSigmaTotComp,dSigmadlambda1=DSigmadlambda1,dSigmadlambda2=DSigmadlambda2,
                            d2Sigmadlambda12=D2SigmaTotdlambda12,Psi=PsiM,Uvec=UTot,Vvec=VTCorrTot,bfuncv=bfuncv,WR=WRTot)
                        
                        d2DDdlambda12<-d2DDdranmat2(X=X,Z=ZOriginal,d2Wdmu2=d2Wdmu2,dWdmu=dWdmu,Wvec=Wvec,dvhatdlambda1=dvhatdlambda1,
                            dvhatdlambda2=dvhatdlambda2,d2vhatdlambda12=d2vhatdlambda12,invSigmaMat=invSigmaTotComp,
                            dSigmadlambda1=DSigmadlambda1,dSigmadlambda2=DSigmadlambda2,d2Sigmadlambda12=D2Sigmadlambda12,WR=WRTot,dWRdu=dWRdUTot,d2WRdu2=d2WRdU2Tot)
                        if (Dindex==1 & Dindex2==1) {
                            d2Vcur12<-d2vhatdlambda12
                            d2Sigm12<-D2Sigmadlambda12
                            d2DDd12<-d2DDdlambda12
                        }                        
                        HessianRVC[firstindex,secondindex]<-d2hdlambda2-0.5*sum(diag(solveDD%*%d2DDdlambda12))+0.5*sum(diag(solveDD%*%dDDdlambda1%*%solveDD%*%dDDdlambda2))
                        
                                      
 
                        # This makes hessian with respect to gamma #
                        HessianRVC[firstindex,secondindex]<-HessianRVC[firstindex,secondindex]*LambdaCur1*LambdaCur2
                    }
                }
            }
        }
    }
    
    # All dispersion parameters must be evaluated jointly #
    # So far we have: HessianODEst - for gammagamma derivatives #
    #                 HessCorr / HessCorrZF - hessian for correlations and fisher z #
    
    
    # Still there could be the same residual variance over the two models that is not implemented yet !!!!!!!!!!!!!!!!!!!!! #
    
    
    
    # Further extension to Truncated Poisson #
    # and another extension to commmon betas etc #
    # further extend to Purahmadi trick #
    # Standard errors and diagnostics #
    # Speed up the algorithm #

    HELPVALUES<-list(X=X,Z=ZOriginal,dWdmu=dWdmu,Wvec=Wvec,dvhatdlambda=dvhatdlambda1,invSigmaMat=invSigmaTotComp,dSigmadlambda=DSigmadlambda1,
                        WR=WRTot,dWRdu=dWRdUTot,y=Y,mu=mu,Phi=Phi,Psi=PsiM,Uvec=UTot,Vvec=VTCorrTot,bfuncv=bfuncv,solveDD=solveDD,VT=VT,
                        dDDd=dDDdlambda1,d2Sigmadlambda12=D2Sigmadlambda12,dVcur1=dVcur1,dVcur2=dVcur2,d2Vcur12=d2Vcur12,dSigm1=dSigm1,dSigm2=dSigm2,
                        d2Sigm12=d2Sigm12,d2DDd12=d2DDd12,dDDd1=dDDd1,dDDd2=dDDd2)
}
if (StandardErrors==FALSE){
    StdErrCorr <- FALSE
    StdErrCorrZF <-FALSE
    StdErrVs <- FALSE
    StdErrBeta <- FALSE
    StdErrODEst <- FALSE
    StdErrDRgamma <- FALSE
}
    OUT<-list(Beta=Beta,V=VT,DRgamma=DRgamma,DYgamma=DYgamma,Correlations=unlist(Correlations),StdErrCorr=StdErrCorr,StdErrCorrZF=StdErrCorrZF,
            StdErrVs=StdErrVs,StdErrBeta=StdErrBeta,StdErrODEst=StdErrODEst,StdErrDRgamma=StdErrDRgamma,
            M2h=M2h,M2pvh=M2pvh,M2pbvh=M2pbvh,CAIC=CAIC)#,HessianOD=HessianODEst,HessianRVC=HessianRVC,FDER=FDER,HessianCorr=HessCorrZF
    if (nrandcor>0) OUT<-c(OUT)
    return(OUT)
}
