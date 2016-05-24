hglmfit_corr <-
function(formulaMain,DataMain,Offset=NULL,RespDist="gaussian",RespLink="identity",
RandDist="gaussian",mord=0,dord=1,spatial=NULL,Neighbor=NULL,Maxiter=200,Iter_mean=5,convergence=10^(-4),
Init_lam=0.25,Init_rho=0.174,contrasts=NULL){
    mc <- match.call()
    fr <- HGLMFrames(mc, formulaMain, contrasts)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formulaMain, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(fr$Y, length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    z<-zz<-z[[1]]
    if (is.null(spatial)) spatial="IND"
    if (spatial=="IAR" && !is.null(Neighbor)) {
          nn<-nrow(Neighbor)
          no<-matrix(0,nn,nn)
          for (i in 1:nn) no[i,i]<-sum(Neighbor[i,])
          sp1<-(no-Neighbor)
          index4<-nn-1
          s_ <-svd(sp1,nu=index4,nv=index4)
          uuu<-matrix(0,nn,index4)
          ddd<-matrix(0,index4,index4)
          for ( i in 1:index4) ddd[i,i]<-sqrt(1/s_$d[i])
          for (i in 1:nn){
             for (j in 1:index4){
             uuu[i,j]<-s_$u[i,j]
          }
          }
          LLL<-uuu %*% ddd
          z<-z%*%LLL
          for (i in 1:nrand) q[i]<-index4
###          print(LLL)
    }
    if (spatial=="MRF_Fix" && !is.null(Neighbor)) {
          rho<-Init_rho
          nn<-nrow(Neighbor)
          no<-matrix(0,nn,nn)
          sp1<-diag(rep(1,nn))
          sp1<-sp1-rho*Neighbor
          index4<-nn
          s_ <-svd(sp1,nu=index4,nv=index4)
          uuu<-matrix(0,nn,index4)
          ddd<-matrix(0,index4,index4)
          for ( i in 1:index4) ddd[i,i]<-sqrt(1/s_$d[i])
          for (i in 1:nn){
             for (j in 1:index4){
             uuu[i,j]<-s_$u[i,j]
          }
          }
          LLL<-uuu %*% ddd
          z<-z%*%LLL
    }
##############################################################
######### initial values : GLM estimates #####################
##############################################################
    dord<-1
    phi <- rep(1,n)
    beta_h<-matrix(0,p,1)
    qcum <- cumsum(c(0, q))
    v_h<-matrix(0,qcum[nrand+1],1)
    if (RandDist=="gaussian") u_h <- v_h
    if (RandDist=="gamma") u_h <-exp(v_h)
    alpha_h <- rep(0, nrand)
    for (i in 1:nrand) alpha_h[i] <- Init_lam
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),offset=Offset)
    beta_h[1:p,1]<-c(resglm$coefficients)[1:p]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<- Offset
    if (spatial=="MRF" || spatial=="MRF_Fix") rho<-Init_rho
    else rho<-0
convergence1<-1
max_iter<-1
while (convergence1>convergence && max_iter<=Maxiter ) {
for(k in 1:Iter_mean) {
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    oq<-matrix(1,qcum[nrand+1],1)
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (spatial=="MRF" || spatial=="MRF_Fix") {
           pW2<-(I-rho*Neighbor)
           W2<-pW2/as.vector(lambda)
    } else {
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }

##############################################################
############# random effect  #################################
##############################################################
    c_v_h<-1.0
    iter_v<-1
    Sig<- z %*% solve(W2) %*% t(z) +solve(W1)
    invSig<-solve(Sig)
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%(y-mu)-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-(1+1/lambda)+exp(-v_h)/lambda
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    v_h_old<-v_h
    v_h<-v_h-solve(d2hdv2)%*%dhdv
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
##    }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    if (mord==0) a<-matrix(0,n,1)
    if (mord==1) {
    T<-t(cbind(t(z),I))
##    Null<-matrix(0,n,n)
##    W<-matrix(0,(2*n),(2*n))
##    W[c(1:n),]<-cbind(W1,Null)
##    W[c((n+1):(2*n)),]<-cbind(Null,W2)   
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
    P<-T%*%solve(t(T)%*%W%*%T)%*%t(T)%*%W
    K1<--z%*%solve(t(T)%*%W%*%T)%*%t(z)
    K2<--solve(t(T)%*%W%*%T)
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    for (i in 1:n){
        d1[i]<-P[i,i]*detadmu[i]
        d2[i]<-0
        for (qq in 1:n){
            d2[i]<-d2[i]+P[qq,qq]*K1[qq,i]
        }
        if (RandDist=="gaussian") d3[i]<-0
    }
    d<-d1+d2+d3
    s<-d*dmudeta/2
    a<-(solve(W1)+z%*%solve(W2)%*%t(z))%*%W1%*%(s*detadmu)
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    beta_h<-solve(t(x)%*%invSig%*%x)%*%(t(x)%*%invSig%*%(z1-a))
    se_beta<-sqrt(diag(solve(t(x)%*%invSig%*%x)))
############################################################## 
} 
###############################################################
############# dispersion parameters ###########################
###############################################################
    v<-v_h
    Q<-invSig-invSig%*%x%*%solve(t(x)%*%invSig%*%x)%*%t(x)%*%invSig
    lam<-alpha_h[1]
#### 1: lam(variance component) , 2: rho
     if (spatial=="MRF" || spatial=="MRF_Fix") {
           pW2<-(I-rho*Neighbor)
           W2<-pW2/as.vector(lambda)
    } else {
        rho<-0
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }
    dREMLdlam<-c(0,0)
    dREML1dlam<-c(0,0)
    dREML2dlam<-c(0,0)
    d2REMLd2lam<-matrix(0,2,2)

    dW2dlam<--W2/lam
    if (spatial=="MRF" && !is.null(Neighbor)) dW2drho<-(-Neighbor)/lam
    else dW2drho<-0

    dSig1dlam<-solve(pW2)
    if (spatial=="MRF" && !is.null(Neighbor)) dSig1drho<--solve(W2)%*%dW2drho%*%solve(W2)
    else dSig1drho<-0

    dvdlam<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2dlam)%*%v
    if (spatial=="MRF" && !is.null(Neighbor)) dvdrho<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho)%*%v
    else dvdrho<-0

    if (RespDist=="gaussian") kkk<-0*mu
    if (RespDist=="poisson") kkk<-mu
    if (RespDist=="binomial") kkk<-(1-2*mu)*dmudeta
    if (RespDist=="gamma") kkk<-0*mu
    dW1dlam<-diag(as.vector(kkk*(z%*%dvdlam)))
    if (spatial=="MRF" && !is.null(Neighbor)) dW1drho<-diag(as.vector(kkk*(z%*%dvdrho)))
    else dW1drho<-0

    dSig2dlam<--solve(W1)%*%dW1dlam%*%solve(W1)
    if (spatial=="MRF" && !is.null(Neighbor)) dSig2drho<--solve(W1)%*%dW1drho%*%solve(W1)
    else dSig2drho<-0

    dSigdlam<-z%*%dSig1dlam%*%t(z)+dSig2dlam
    if (spatial=="MRF" && !is.null(Neighbor)) dSigdrho<-z%*%dSig1drho%*%t(z)+dSig2drho
    else  dSigdrho<-0

    dterm1dlam<--t(v)%*%(dW2dlam)%*%v/2
    if (spatial=="MRF" && !is.null(Neighbor)) dterm1drho<--t(v)%*%(dW2drho)%*%v/2
    else dterm1drho<-0

    if (RespDist=="poisson") dW1dv<-W1
    if (RespDist=="gaussian") dW1dv<-diag(as.vector(kkk))
    if (RespDist=="binomial") dW1dv<-diag(as.vector(kkk))
    if (RespDist=="gamma") dW1dv<-diag(as.vector(kkk))
    dterm2dv<-y-mu-z%*%W2%*%v-1/2*(1/diag(W1))*diag(dW1dv)

    dREMLdlam[1]<--t(v)%*%dW2dlam%*%v/2-t(dvdlam)%*%W2%*%v-0.5*sum(diag(Q%*%dSigdlam))+t(dvdlam)%*%t(z)%*%W1%*%((y-mu)*detadmu)-0.5*sum(diag(solve(W1)%*%dW1dlam))
    if (spatial=="MRF" && !is.null(Neighbor)) dREMLdlam[2]<-t(dvdrho)%*%t(z)%*%W1%*%((y-mu)*detadmu)-0.5*sum(diag(solve(W1)%*%dW1drho))-t(dvdrho)%*%W2%*%v-t(v)%*%dW2drho%*%v/2 -0.5*sum(diag(Q%*%dSigdrho))

    d2W2dlam2<-2*W2/lam^2
    d2W2drho2<-matrix(0,q,q)
    if (spatial=="MRF" && !is.null(Neighbor)) d2W2dlamrho<-Neighbor/lam
    d2vdlam2<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2dlam%*%dvdlam-d2W2dlam2%*%v)
    if (spatial=="MRF" && !is.null(Neighbor)) d2vdrho2<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho%*%dvdrho-d2W2drho2%*%v)

    H<-t(z)%*%W1%*%z+W2
    dHdlam<-dW2dlam
    if (spatial=="MRF" && !is.null(Neighbor)) dHdrho<-dW2drho
    d2Hdlam2<-d2W2dlam2
    if (spatial=="MRF" && !is.null(Neighbor)) d2Hdrho2<-d2W2drho2
    if (spatial=="MRF" && !is.null(Neighbor)) d2Hdlamrho<-d2W2dlamrho

    d2REMLd2lam[1,1]<--0.5*t(v)%*%d2W2dlam2%*%v-0.5*sum(diag(solve(W2)%*%dW2dlam%*%solve(W2)%*%dW2dlam))+0.5*sum(diag(solve(W2)%*%d2W2dlam2))+0.5*sum(diag(solve(H)%*%dHdlam%*%solve(H)%*%dHdlam))-0.5*sum(diag(solve(H)%*%d2Hdlam2)) 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[1,2]<--0.5*t(v)%*%d2W2dlamrho%*%v-0.5*sum(diag(solve(W2)%*%dW2dlam%*%solve(W2)%*%dW2drho))+0.5*sum(diag(solve(W2)%*%d2W2dlamrho))+0.5*sum(diag(solve(H)%*%dHdlam%*%solve(H)%*%dHdrho))-0.5*sum(diag(solve(H)%*%d2Hdlamrho)) 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[2,1]<-d2REMLd2lam[1,2] 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[2,2]<-0-0.5*sum(diag(solve(W2)%*%dW2drho%*%solve(W2)%*%dW2drho))+0.5*sum(diag(solve(W2)%*%d2W2drho2))+0.5*sum(diag(solve(H)%*%dHdrho%*%solve(H)%*%dHdrho))-0.5*sum(diag(solve(H)%*%d2Hdrho2))  

    clam<-c(lam,rho)
    old_clam<-clam

    if (spatial=="MRF" && !is.null(Neighbor)) {
        clam<-clam-solve(d2REMLd2lam)%*%dREMLdlam
##        print(clam)
##        print(dREMLdlam)
##        print(d2REMLd2lam)
        if (clam[2]>1) clam[2]<- clam[2]
        if (clam[2]< -1) clam[2]<- clam[2]
    }
    else clam[1]<-clam[1]-dREMLdlam[1]/d2REMLd2lam[1,1]
##    print(clam[2])
    convergence1<-sum(abs(clam-old_clam))
    lam<-clam[1]
    rho<-clam[2]
    if (spatial=="MRF_Fix") rho<-Init_rho
    alpha_h[1]<-lam
    max_iter<-max_iter+1
    print_i<-max_iter
    print_err<-convergence1
    names(print_i) <- "iteration : "
##    print(print_i)
    names(print_err) <- "convergence : "
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (spatial=="MRF" || spatial=="MRF_Fix") {
        pW2<-(I-rho*Neighbor)
        W2<-pW2/as.vector(lambda)
    } else {
        rho<-0
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }
    if (RespDist=="gaussian") {
        temp5<- dmudeta^2 /(Vmu)
        dW1dphi<-diag(-as.vector(temp5)/phi^2)
        dvdphi<-solve(t(z)%*%W1%*%z+W2)%*%(-t(z)%*%dW1dphi%*%(y-mu))
###        dvdphi<-0*dvdphi
        dSigdphi<-diag(as.vector(temp5))
        dW2dphi<-0*dW2dlam
        d2W2dphi2<-0*d2W2dlam2
        H<-t(z)%*%W1%*%z+W2
        HX<-t(x)%*%W1%*%x
        HXZ<-t(x)%*%W1%*%z
        dHdphi<-t(z)%*%dW1dphi%*%z
        dHXdphi<-t(x)%*%dW1dphi%*%x
        dHXZdphi<-t(x)%*%dW1dphi%*%z
        Hp<-rbind(cbind(HX,HXZ),cbind(t(HXZ),H))
        dHpdphi<-rbind(cbind(dHXdphi,dHXZdphi),cbind(t(dHXZdphi),dHdphi))
        d2W1dphi2<-diag(2*as.vector(temp5)/phi^3)
        d2Hdphi2<-t(z)%*%d2W1dphi2%*%z
        d2HXdphi2<-t(x)%*%d2W1dphi2%*%x
        d2HXZdphi2<-t(x)%*%d2W1dphi2%*%z
        d2Hpdphi2<-rbind(cbind(d2HXdphi2,d2HXZdphi2),cbind(t(d2HXZdphi2),d2Hdphi2))
        dREMLdphi<--t(v)%*%dW2dphi%*%v/2-t(dvdphi)%*%W2%*%v-0.5*sum(diag(Q%*%dSigdphi))+t(dvdphi)%*%t(z)%*%W1%*%((y-mu)*detadmu)+sum(0.5*(y-mu)^2/phi^2)
###            -0.5*sum(diag(solve(W1)%*%dW1dphi))
###        d2REMLd2phi<--0.5*t(v)%*%d2W2dphi2%*%v-0.5*sum(diag(solve(W2)%*%dW2dphi%*%solve(W2)%*%dW2dphi))+0.5*sum(diag(solve(W2)%*%d2W2dphi2))-sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2)) 
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)-sum((p+q)*0.5/(n*phi^2))
###        print(dREMLdphi)
###        print(d2REMLd2phi)
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2))+0.5*sum(diag(solve(HX)%*%dHXdphi%*%solve(HX)%*%dHXdphi))-0.5*sum(diag(solve(HX)%*%d2HXdphi2)) 
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi)) 
        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi))-0.5*sum(diag(solve(Hp)%*%d2Hpdphi2))
        oldphi1<-phi[1]
        phi1<-phi[1]+dREMLdphi/abs(d2REMLd2phi)
        phi <- rep(phi1,n)
        convergence2<-sum(abs(phi[1]-oldphi1))
        convergence1<-convergence1+convergence2
        print_err<-print_err+convergence2
###        print(phi[1])
    }
    if (RespDist=="gamma") {
        temp5<- dmudeta^2 /(Vmu)
        dW1dphi<-diag(-as.vector(temp5)/phi^2)
        dvdphi<-solve(t(z)%*%W1%*%z+W2)%*%(-t(z)%*%dW1dphi%*%(y-mu))
###        dvdphi<-0*dvdphi
        dSigdphi<-diag(as.vector(temp5))
        dW2dphi<-0*dW2dlam
        d2W2dphi2<-0*d2W2dlam2
        H<-t(z)%*%W1%*%z+W2
        HX<-t(x)%*%W1%*%x
        HXZ<-t(x)%*%W1%*%z
        dHdphi<-t(z)%*%dW1dphi%*%z
        dHXdphi<-t(x)%*%dW1dphi%*%x
        dHXZdphi<-t(x)%*%dW1dphi%*%z
        Hp<-rbind(cbind(HX,HXZ),cbind(t(HXZ),H))
        dHpdphi<-rbind(cbind(dHXdphi,dHXZdphi),cbind(t(dHXZdphi),dHdphi))
        d2W1dphi2<-diag(2*as.vector(temp5)/phi^3)
        d2Hdphi2<-t(z)%*%d2W1dphi2%*%z
        d2HXdphi2<-t(x)%*%d2W1dphi2%*%x
        d2HXZdphi2<-t(x)%*%d2W1dphi2%*%z
        d2Hpdphi2<-rbind(cbind(d2HXdphi2,d2HXZdphi2),cbind(t(d2HXZdphi2),d2Hdphi2))
###        d2HXdphi2<-0*t(x)%*%dW1dphi%*%x
###        dREMLdphi<-sum(-log(y)/phi^2+y/(phi^2*mu)+log(phi)/phi^2-1/phi^2+log(mu)/phi^2+digamma(1/phi)/phi^2)-0.5*sum(diag(solve(H)%*%dHdphi))-0.5*sum(diag(solve(HX)%*%dHXdphi))
        dREMLdphi<-sum(-log(y)/phi^2+y/(phi^2*mu)+log(phi)/phi^2-1/phi^2+log(mu)/phi^2+digamma(1/phi)/phi^2)-0.5*sum(diag(solve(Hp)%*%dHpdphi))
###        d2REMLd2phi<-sum(2*log(y)/phi^3-2*y/(phi^3*mu)-2*log(phi)/phi^3+3/phi^3-log(mu)/phi^3-2*digamma(1/phi)/phi^3-trigamma(1/phi)/phi^4)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2))+0.5*sum(diag(solve(HX)%*%dHXdphi%*%solve(HX)%*%dHXdphi))
        d2REMLd2phi<-sum(2*log(y)/phi^3-2*y/(phi^3*mu)-2*log(phi)/phi^3+3/phi^3-2*log(mu)/phi^3-2*digamma(1/phi)/phi^3-trigamma(1/phi)/phi^4)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi))-0.5*sum(diag(solve(Hp)%*%d2Hpdphi2))
        oldphi1<-phi[1]
        phi1<-phi[1]+dREMLdphi/abs(d2REMLd2phi)
        phi <- rep(phi1,n)
        convergence2<-sum(abs(phi[1]-oldphi1))
        convergence1<-convergence1+convergence2
        print_err<-print_err+convergence2
    }
##    print(print_err)
}
###############################################################
############# se for dispersion estimates######################
###############################################################
    X<-x
    p<-ncol(X)
    O1<-matrix(0,p,p)
    O2<-matrix(0,p,qcum[nrand+1])
    if (spatial=="MRF" && !is.null(Neighbor)) infoterm<-matrix(0,2,2)
    else infoterm<-matrix(0,2,2) 
    d2hlikedlam2<-n/(2*lam^2)-t(v)%*%pW2%*%v/lam^3
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%z%*%I)),cbind((t(I)%*%t(z)%*%W1%*%X),H))
    dAdlam<-rbind(cbind(O1,O2),cbind(t(O2),dHdlam))
    d2Adlam2<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdlam2))
    d2hlikedlamdv<-pW2%*%v/lam^2
    dAdv_dvdlam<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%z%*%I)),
                 cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%z%*%I)))
    dAdv_d2vdlam2<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%z%*%I)),
                 cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%z%*%I)))
    if (spatial=="MRF" && !is.null(Neighbor))  {
        d2hlikedrho2<--1/2*sum(diag(solve(pW2)%*%Neighbor%*%solve(pW2)%*%Neighbor))
        d2hlikedrhodlam<--t(v)%*%Neighbor%*%v/(2*lam^2)
        d2vdrhodlam<--solve(H)%*%dHdrho%*%solve(H)%*%(pW2%*%v)/lam^2+solve(H)%*%(-Neighbor%*%v)/lam^2
        dAdrho<-rbind(cbind(O1,O2),cbind(t(O2),dHdrho))
        d2Adrho2<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdrho2))
        d2Adrhodlam<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdlamrho))
        d2hlikedrhodv<--Neighbor%*%v/lam
        dAdv_dvdrho<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%z%*%I)))
        dAdv_d2vdrho2<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%z%*%I)))
        dAdv_d2vdrhodlam<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%z%*%I)))
   }

    tinfoterm1<-d2hlikedlam2+1/2*sum(diag(solve(A)%*%dAdlam%*%solve(A)%*%dAdlam))-1/2*sum(diag(solve(A)%*%d2Adlam2))
    tinfoterm2<-sum(as.vector(d2hlikedlamdv*dvdlam))+1/2*sum(diag(solve(A)%*%(dAdv_dvdlam)%*%solve(A)%*%dAdlam))
    tinfoterm3<-tinfoterm2-1/2*sum(diag(solve(A)%*%dAdv_d2vdlam2))
    tinfoterm4<-t(dvdlam)%*%(-H)%*%dvdlam+1/2*sum(diag(solve(A)%*%dAdv_dvdlam%*%solve(A)%*%dAdv_dvdlam))

    infoterm[1,1]<-tinfoterm1+tinfoterm2+tinfoterm3+tinfoterm4
    if (spatial=="MRF" && !is.null(Neighbor))  {    
         rinfoterm1<-d2hlikedrho2+1/2*sum(diag(solve(A)%*%dAdrho%*%solve(A)%*%dAdrho))-1/2*sum(diag(solve(A)%*%d2Adrho2))
         rinfoterm2<-sum(as.vector(d2hlikedrhodv*dvdrho))+1/2*sum(diag(solve(A)%*%(dAdv_dvdrho)%*%solve(A)%*%dAdrho))
         rinfoterm3<-rinfoterm2-1/2*sum(diag(solve(A)%*%dAdv_d2vdrho2))
         rinfoterm4<-t(dvdrho)%*%(-H)%*%dvdrho+1/2*sum(diag(solve(A)%*%dAdv_dvdrho%*%solve(A)%*%dAdv_dvdrho))
         infoterm[2,2]<-rinfoterm1+rinfoterm2+rinfoterm3+rinfoterm4
    trinfoterm1<-d2hlikedrhodlam+1/2*sum(diag(solve(A)%*%dAdrho%*%solve(A)%*%dAdlam))-1/2*sum(diag(solve(A)%*%d2Adrhodlam))
    trinfoterm2<-sum(as.vector(d2hlikedlamdv*dvdrho))+1/2*sum(diag(solve(A)%*%(dAdv_dvdrho)%*%solve(A)%*%dAdlam))
    trinfoterm3_1<-sum(as.vector(d2hlikedrhodv*dvdlam))+1/2*sum(diag(solve(A)%*%(dAdv_dvdlam)%*%solve(A)%*%dAdrho))
    trinfoterm3<-trinfoterm3_1-1/2*sum(diag(solve(A)%*%dAdv_d2vdrho2))
    trinfoterm4<-t(dvdrho)%*%(-H)%*%dvdlam+1/2*sum(diag(solve(A)%*%dAdv_dvdrho%*%solve(A)%*%dAdv_dvdlam))
    infoterm[1,2]<-infoterm[2,1]<-trinfoterm1+trinfoterm2+trinfoterm3+trinfoterm4
    }
    clam_se<-matrix(0,2,1)
    if (spatial=="MRF" && !is.null(Neighbor)) {
         temp4<-sqrt(abs(diag(solve(-infoterm))))
         for (i in 1:2) clam_se[i,1]<-temp4[i]
    }
    else clam_se[1,1]<-sqrt(abs(-1/infoterm[1,1]))
###############################################################
############# likelihood estimates ############################
###############################################################
    pi<-3.14159265359
    d2hdv2<--t(z)%*%W1%*%z-W2
    H<-t(z)%*%W1%*%z+W2
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%z%*%I)),cbind((t(I)%*%t(z)%*%W1%*%X),H))
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/phi-0.5*log(2*pi*phi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-log(factorial(y)))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu)+(1-y)*log(1-mu))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/phi-log(y)-y/(phi*mu)-log(phi)/phi-log(mu)/phi-lgamma(1/phi))
    hlikeli1<--2*hlikeli
    hlikeli<-hlikeli-0.5*t(v_h)%*%W2%*%v_h-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))
    pvh<-hlikeli-0.5*log(abs(det(-d2hdv2/(2*pi))))
    pbvh<-hlikeli-0.5*log(abs(det(A/(2*pi))))
    m2h<--2*hlikeli
    m2pvh<--2*pvh
    m2pbvh<--2*pbvh
###############################################################
############# print estimates ###########################
###############################################################
    z_beta<-beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff<-cbind(beta_h,se_beta,z_beta)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
##    print("Estimates from the mean model")    
##    print(beta_coeff,4)
    if (RespDist=="gaussian") {    
##        print("Estimates from the dispersion model for Phi")    
        se_phi<-sqrt(-1/d2REMLd2phi)
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
##        print(phi_coeff,4)
    }
    if (RespDist=="gamma") {    
##        print("Estimates from the dispersion model for Phi")    
        se_phi<-sqrt(-1/d2REMLd2phi)
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
##        print(phi_coeff,4)
    }
##    if (spatial=="IAR" && !is.null(Neighbor)) print("Estimates from the dispersion model for Lambda in the IAR model")
##    else print("Estimates from the dispersion model for Lambda")
    se_lam<-clam_se[1,1]
    z_lam<-lam/se_lam
    lam_coeff<-cbind(lam,se_lam)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
##     print(lam_coeff,4)
    if (spatial=="MRF" && !is.null(Neighbor)) {
##        print("Estimates for rho in the MRF model")
        se_rho<-clam_se[2,1]
        z_rho<-rho/se_rho
        rho_coeff<-cbind(rho,se_rho)
        colnames(rho_coeff) <- c("Estimate", "Std. Error")
        rownames(rho_coeff) <- "rho"
##        print(rho_coeff,4)
    } else se_rho<-0.0001
###############################################################
############# Likelihoods         ###########################
###############################################################
    if (dord<=1) like_value<-cbind(m2h,m2pvh,m2pbvh)
    if (dord<=1) colnames(like_value) <- c("-2*h","-2*p_v(h)","-2p_b,v(h)")
##    print(like_value)
    res<-list(namesX,beta_h,se_beta,lam,rho,clam_se,v_h,like_value,hlikeli1,mu,W2,se_lam,se_rho,A)
    return(res)
}
