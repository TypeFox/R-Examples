############## Two-Level Normal independent sampling estimation  ##############
##############       SUPPORTS MULTIVARIATE OUTCOMES
############## Copyright 2000, Phil Everson, Swarthmore College. ##############

######################################################################
## Minor changes for R port Copyright (C) 2004-2005, Roger D. Peng <rpeng@jhsph.edu>
######################################################################
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###############################################################################


tlnise<-function(Y,V,w=NA,V0=NA,prior=NA,N=1000,seed=NULL,
                 Tol=1e-6,maxiter=1000,intercept=TRUE,labelY=NA,labelYj=NA,
                 labelw=NA,digits=4,brief=1,prnt=TRUE){
    ##
    ## This program is free and may be redistributed as long as
    ## the copyright is preserved. 
    ##
    ## This program should be cited if used in published research.
    ##
    ## The author will appreciate notification of any comments or
    ## errors by emailing peverso1@swarthmore.edu
    ##
    ## Reference: 
    ## Everson and Morris (2000). J. R. Statist. Soc. B, 62 prt 2, pp.399-412.
    ##
    ## 			!!!IMPORTANT!!! 
    ## you must change the path in the following line to specify the directory in
    ## which you stored the Fortran object file tlnise.o.
    ## dyn.load("tlnisemv1.o")
    ## dyn.load("tlnisemv1.so")
    ## dyn.open("S.so")
    ##          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ -> change to agree with your path
    ##
    ## TLNise: Two-Level Normal independent sampling estimation 
    ##
    ##  Provides inference for 2-level Normal hierarchical models 
    ##  with univariate (p = 1) or multivariate (p > 1) outcomes, and
    ##  with known level-1 variances V:
    ##
    ##  Level-1:  Yj|thetaj, Vj ~ Np(thetaj, Vj), j = 1,..,J
    ##
    ##  Level-2:  thetaj|A, gamma ~ Np(Wjgamma,A)
    ##
    ##  Calls Fortran subroutines in tlnise.f (object file tlnise.o).
    ##
    ## INPUTS:
    ##        Y: Jxp (or pxJ) matrix of p-dimensional Normal outcomes
    ##        V: pxpxJ array of pxp Level-1 covariances (assumed known)
    ##        w: Jxq (or qxJ) covariate matrix (adds column of 1's if 
    ##	     not included and intercept = T)
    ##       V0: "typical" Vj (default is average of Vj's)
    ##    prior: Prior parameter: p(B0) = |B0|^{(prior-p-1)/2}
    ##      prior = -(p+1): uniform on level-2 covariance matrix A (default) 
    ##      prior =      0: Jeffreys' prior
    ##      prior =  (p+1): uniform on shrinkage matrix B0
    ##        N: number of Constrained Wishart draws for inference
    ##      Tol: tolerance for determining modal convergence
    ##  maxiter: maximum number of EM iterations for finding mode
    ##intercept: if True, an intercept term is included in the regression
    ##   labelw: optional names vector for covariates
    ##   labelY: optional names vector for the J observations
    ##  labelYj: optional names vector for the p elements of Yj
    ##   digits: number of significant digits for reporting results
    ##    brief: level of output, from 0 (minimum) to 2 (maximum)
    ##     prnt: Controls printing during execution
    ##
    ## OUTPUTS:
    ##  *** brief=0 ***
    ##    gamma:  matrix of posterior mean and sd estimates of Gamma,
    ##             and their ratios.
    ##    theta:  pxJ matrix of posterior mean estimates for thetaj's.
    ##  SDtheta:  pxJ matrix of posterior SD estimates for thetaj's.
    ##        A:  pxp estimated posterior mean of variance matrix A.
    ##      rtA:  p-vector of between group SD estimates.
    ##
    ##  *** brief=1 *** 
    ##   brief=0 outputs, plus
    ##  Dgamma: rxr estimated posterior covariance matrix for Gamma.
    ##  Vtheta: pxpxJ array of estimated covariances for thetaj's.
    ##      B0: pxpxN array of simulated B0 values.
    ##      lr: N-vector of log density ratios for each B0 value.
    ##  
    ##  *** brief=2 *** 
    ##   brief=1 outputs, plus
    ##      lf: N-vector of log f(B0|Y) evaluated at each B0.
    ##     lf0: N-vector of log f0(B0|Y) evaluated at each B0
    ##           (f0 is the CWish envelope density for f).
    ##      df: degrees of freedom for f0.
    ##   Sigma: scale matrix for f0.
    ##    nvec: number of matrices begun, diagonal and off-diagonal
    ##           elements simulated to get N CWish matrices.
    ##    nrej: number of rejections that occurred at each step 1,..,p.
######################
    ## Check inputs:
    out.chk<-checkcon(Y,V,w,intercept,prior,prnt)
    ## sets Y: pxJ, w: qxJ, V: pxpxJ
    Y<-out.chk$Y
    V<-out.chk$V
    w<-out.chk$w
    J<-out.chk$J
    p<-out.chk$p
    q<-out.chk$q
    r<-p*q
    prior<-out.chk$prior
    if(prnt)
        print(paste("******** Prior Parameter =",prior),quote=FALSE)
    if(missing(V0))
            V0 <- rowMeans(V, dims = 2)
    if(max(abs(V0-diag(p))) < Tol){
        Ys<-Y
        Vs<-V
        rtV0<-diag(p)
    }
    else{
        ## Rotate Y and V by rtV0^{-1}:
        newvars<-standard.f(Y,V,V0)
        Ys<-newvars$Y
        Vs<-newvars$V
        rtV0<-newvars$rtVo
        ## Ys = rtV0^-1%*%Y, Vsj = rtV0^(-1)%*%V%*%rtV0^(-1)
    }
    ##
    ## Locate A corresponding to the posterior mode of 
    ##  B0 = (I+A*)^{-1}  (A* = rtV0^-1%*%A%*%rtV0^-1).
    if(prnt){
        cat("\n")
        print("******** Locating Posterior Mode ************ ",quote=FALSE)
        cat("\n")
    }
    Astart<-diag(p)
    out.mode<-postmode.f(Ys,Vs,w,rtV0,prior,Astart,Tol, maxiter,J,p,q,r)
    modeB0<-solve(diag(p)+out.mode$newA)
    lfmode<-out.mode$lf
    ## compute A value on scale of data:
    modeA<-rtV0%*%out.mode$newA%*%rtV0
    ##
    ## Set the degrees of freedom (df) and scale parameter (Sigma) for
    ##  the CWish density f0 so that the mode of f0 is the same as  
    ##  the mode of f(B0|Y):
    df<-J-q+prior
    Sigma<-modeB0/(df-p-1)
    ## mode(f0) = (df-p-1)*Sigma = modeB0.
    ##
    ## Set d as the eigenvalues of Sigma^{-1}, ordered smallest to largest:
    eigS<-eigen(Sigma,symmetric=TRUE)
    d<-1/eigS[[1]]
    ##
    ## Compute the inverse of Sigma, Siginv, and a square root matrix
    ##  for Sigma, rtSig ( rtSig%*%t(rtSig) = Sigma).
    if(p==1){
        Siginv<-matrix(1/Sigma)
        rtSig<-matrix(sqrt(Sigma))
    }
    else{
        Siginv<-eigS[[2]]%*%diag(d)%*%t(eigS[[2]])
        rtSig<-eigS[[2]]%*%diag(sqrt(1/d))
    }
    ## 
    ## Compute the log of the CWish(df,Sigma;d) density at its mode:
    lf0mode<-(df-p-1)*ldet(modeB0)/2 - tr(modeB0%*%Siginv)/2
    ##
    ## Set adj as the difference in the log-densities at the mode. This
    ## will be used to adjust the log-density ratio before exponentiating. 
    adj<-lf0mode-lfmode
    iter<-out.mode[[16]]
    if(prnt){
        if(iter<maxiter){
            print(paste("Converged in",iter,"EM iterations."),quote=FALSE)
        }
        else{
            print(paste("Did not converge in maxiter =",maxiter,"EM steps."),quote=FALSE)
        }
        print("Posterior mode of B0:",quote=FALSE)
        cat("\n")
        print(modeB0,digits)
        cat("\n")
        print(paste(
                    "lf(modeB0) =",signif(lfmode,digits),"; lf0(modeB0) =",
                    signif(lf0mode,digits),"; adj =", signif(adj,digits)),quote=FALSE)
    }
    if(is.null(seed))
            seed <- ceiling(runif(1) * 1e8)
    if(seed > 0)
            seed <- -seed
    if(prnt) {
        cat("\n")
        print("******** Drawing Constrained Wisharts ******** ",quote=FALSE)
        cat("\n")
    }
    ## Generate N Constrained-Wishart(df,I;D) draws:
    pd<-pchisq(d,df)
    if(min(pd)> 1-Tol){
        pd<-rep(1,p)
        d<-rep(0,p)
    }
    outU<-rscwish(N,p,df,d,pd,seed)
    Uvals<-outU$Ua
    nvec<-outU$nvec
    nrej<-outU$nrej
    if(prnt){
        cat("\n")
        print(paste("CWish acceptance rate =",N,"/",nvec[1],"=",
                    signif(N/nvec[1],digits)),quote=FALSE)
    }
    ##
    ## Rotate Ui's to get B0i's:
    B0vals<-mammult(rtSig, Uvals, t(rtSig))
    ##
    if(prnt){
        cat("\n")
        print("******** Processing Draws ************** ",quote=FALSE)
        cat("\n")
    }
    out.lf<-lfB0.f(B0vals,Ys,Vs,w,rtV0,df,Siginv,N,J,p,q,r,prior,adj)
    ##
    lr<-out.lf$lrv
    lf<-out.lf$lfv
    lf0<-out.lf$lf0v
    ## changed 7/31/2000 for brief=2 output.
    avewt<-mean(exp(lr-max(lr)))
    if(prnt){
        print(paste("Average scaled importance weight =",
                    signif(avewt,digits)),quote=FALSE)
        cat("\n")
    }
    meanA<-rtV0%*%(out.lf$meanA)%*%rtV0
    if(p==1){ rtA<-sqrt(meanA)}
    else{ rtA<-sqrt(diag(meanA))}
    ##
    if(prnt){
        print("Posterior mean estimate of A:",quote=FALSE)
        print(meanA, digits)
        cat("\n")
        print("Between-group SD estimate:", quote=FALSE)
        print(rtA,digits)
        cat("\n")
    }
    ## Set theta:
    if(missing(labelY)) labelY<-1:J
    if(missing(labelYj)) labelYj<-1:p
    theta<-rtV0%*%out.lf$thetahat
    Vtheta<-mammult(rtV0,out.lf$Vthetahat,rtV0)
    if(p==1){ 
        SDtheta<-sqrt(c(Vtheta))
        theta<-c(theta)
        names(SDtheta)<-names(theta)<-labelY
    }
    else{
        SDtheta<-0*Y
        for(j in 1:J) SDtheta[,j]<-sqrt(diag(Vtheta[,,j]))
        dimnames(theta)<-dimnames(SDtheta)<-list(labelYj,labelY)
        dimnames(Vtheta)<-list(labelYj,labelYj,labelY)
    }
    ##
    gamma<-out.lf$gamhat
    Dgamma<-out.lf$Dgamhat
    GammaMat<-cbind(gamma,sqrt(diag(Dgamma)))
    GammaMat<-cbind(GammaMat,GammaMat[,1]/GammaMat[,2])
    if(missing(labelw)) labelw<-seq(0,r-1,1)
    dimnames(GammaMat)<-list(labelw, c("est","se","est/se"))
    names(gamma)<-labelw
    dimnames(Dgamma)<-list(labelw,labelw)
    ## Output:
    if(brief==0)
        out<-list(gamma=GammaMat, theta=theta, SDtheta=SDtheta, A=meanA, rtA=rtA)
    if(brief==1)
        out<-list(gamma=GammaMat, theta=theta, SDtheta=SDtheta, A=meanA, rtA=rtA, 
                  Dgamma=Dgamma,Vtheta=Vtheta, B0=B0vals, lr=lr)
    if(brief==2)
        out<-list(gamma=GammaMat, theta=theta, SDtheta=SDtheta, A=meanA, rtA=rtA, 
                  Dgamma=Dgamma,Vtheta=Vtheta, B0=B0vals, lr=lr,lf=lf,
                  lf0=lf0,df=df,Sigma=Sigma,nvec=nvec,nrej=nrej)
    if(brief>2)
        out<-list(out.mode, outU=outU, out.lf=out.lf)
    out}


postmode.f<-function(Y,V,w,rtV0,prior,Astart,EPS,MAXIT,K,p,q,r){
    ## Calls Fortran subroutine postmode to locate the A value
    ## corresponding to the posterior mode of f(B0|Y).
    ## Y: pxJ, V: pxpxJ, w: qxJ
    indxr<-gamma<-rep(0,r)
    indxp<-rep(0,p)
    indxpk<-0*Y
    H<-newA<-0*Astart
    Dgam<-matrix(0,ncol=r,nrow=r)
    lf<-llik<-0
    Wi<-matrix(0,ncol=p,nrow=r)
    out<-.Fortran("postmode",
                  Astart=Astart,as.double(Y), as.double(w),as.double(V),
                  as.double(rtV0),as.double(prior), ## 6
                  as.integer(p), as.integer(q),as.integer(r), as.integer(K),
                  as.double(EPS), as.integer(MAXIT),    ## 12
                  newA=newA,
                  lf=lf,
                  llik=llik,
                  as.integer(1), ## 16 iter
                  as.double(gamma),as.double(indxpk),
                  as.double(Dgam),as.double(Dgam),
                  as.double(0*V), 
                  as.double(indxp),  ## Yi
                  as.double(Wi),   ## Wi  23
                  as.double(t(Wi)), 
                  as.double(H),  ## Di   
                  as.double(H),  ## Si
                  H=H,  ## H
                  as.double(H),  ## Hi
                  as.integer(indxp), 
                  as.integer(indxp), ## indxp2 
                  as.integer(indxpk),   ## 30
                  as.integer(indxr), PACKAGE = "tlnise")
    out
}


rscwish<-function(N,p,df,d,pd,seed){
    ## Calls Fortran subroutine "rscwish" to generate N pxp
    ## standard constrained Wishart matrices with df degrees of
    ## freedom and a p-vector of constraints d. pd is the p-vector
    ## of the Chi-square(df) CDF evaluated at d. If all elements of 
    ## d are 0 and all elements of pd are 1, then rscwish returns N 
    ## unconstrained Wisharts.
    ##
    Ua<-array(0.0,c(p,p,N))
    Tmat<-matrix(0.0,p,p)
    Z<-rep(0,p)
    chitab<-c(qchisq(seq(.001,.999,.001), df))
    ## chitab is a table of 999 Chi-square(df) quantiles
    ##  to pass to Fortran. 
    out<-.Fortran("rscwish",
                  Ua=Ua,
                  as.integer(N), 
                  as.integer(p), 
                  as.double(df),
                  as.double(d), 
                  as.double(pd), 
                  as.double(chitab),
                  as.integer(seed),
                  as.double(Tmat),     ## Tmat
                  as.double(Tmat),     ## Deltinv
                  as.double(Tmat),     ## H11
                  as.double(Z),        ## H12
                  as.double(Z),        ## Zt
                  as.double(Z),        ## Zts
                  as.double(Z),        ## tempp
                  as.integer(0),       ## mcnt
                  as.integer(0),       ## drvcnt
                  as.integer(0),       ## orvcnt
                  as.integer(Z),        ## rejcnt
                  PACKAGE = "tlnise")
    list(Ua=out[[1]],nvec=c(out[[16]],out[[17]],out[[18]]),nrej=c(out[[19]]))
}



lfB0.f<-function(B0vals, Y, V, w, rtV0, df,Siginv,N,J,p,q,r, 
                 prior=NA,adj=0){
    ## Returns the log posterior density of B0 = (I+A)^{-1},
    ## assuming prior distribution:
    ##   p(B0)dB0 propto |B0|^{(prior-p-1)/2} dB0, 0 < B0 < I.
    ## Assumes data have been rotated by rtV0^-1.
    if(missing(prior))
        ## Uniform prior on B0 -> reports Likelihood of B0.
        prior<-p+1
    lfv<-lf0v<-lrv<-rep(0,N)
    gamma<-rep(0,r)
    Dgam<-matrix(0,ncol=r,nrow=r)
    theta<-0*Y
    Vtheta<-0*V
    Yj<-rep(0,p)
    Wj<-matrix(0,ncol=r,nrow=p)
    A<-matrix(0,ncol=p,nrow=p)
    resid<-0*Y
    sumwt<-0
    out<-.Fortran("lfb0", B0vals=B0vals,as.double(Y),as.double(V), 
                  as.double(w), as.double(rtV0), 
                  as.integer(N),as.integer(p),as.integer(q),as.integer(r),   ## 9
                  as.integer(J),as.double(df),as.double(Siginv), ##12
                  as.double(prior),as.double(adj),lfv=lfv,lf0v=lf0v,lrv=lrv,  ##17
                  gamhat=gamma,Dgamhat=Dgam, thetahat=theta, Vthetahat=Vtheta, 
                  meanA=A,sumwt=sumwt,
                  as.double(gamma),as.double(Dgam), as.double(Dgam), as.double(Yj),  ##27
                  as.double(Vtheta),as.double(resid),as.double(A),as.double(A),  
                  as.double(A),as.double(Dgam), as.double(Vtheta),as.double(Yj), ##35
                  as.double(A),as.double(Wj), as.double(t(Wj)),as.double(A),
                  as.integer(Yj),as.integer(Yj),as.integer(resid),as.integer(gamma), ##43
                  PACKAGE = "tlnise")
    out
}


checkcon<-function(Y,V,w,intercept,prior,prnt){
    if(length(Y)==length(V)){
        ## univariate problem:
        J<-length(Y)
        Y<-matrix(Y,ncol=J,nrow=1)
        V<-array(V,c(1,1,J))
    }
    dmV<-dim(V)
    p<-dmV[1]
    if(p!=dmV[2]) stop("error in dimension of V")
    J<-dmV[3]
    dmY<-dim(Y)
    if(dmY[1]==J){
        if(dmY[2]!=p) stop("error in dimension of Y or V")
        ## Set Y so it is pxJ:
        Y<-t(Y)
    }
    else{
        if(dmY[1]!=p|dmY[2]!=J) stop("error in dimension of Y or V")
    }
    if(is.na(w[1])){
        w<-matrix(rep(1,J),nrow=1)
        q<-1
    }
    else{
        if(length(w)==J)
            w<-matrix(w,nrow=1)
        dmw<-dim(w)
        if(dmw[1]==J){
            q<-dmw[2]
            w<-t(w)
        }
        else{
            if(dmw[2]!=J) stop("error in dimension of w or V")
            q<-dmw[1]
        }
    }
    if(intercept&&max(abs(w[1,]-1))!=0){
        ##  add a row of 1's for an intercept:
        w<-rbind(rep(1,J),w)
        q<-q+1
    }
    ## Returns w: qxJ
    if(J-q < 1) 
        stop("too many covariates for the number of observations")
    if(!is.na(prior)) {
        if(J-q+prior < 1){
            if(prnt){
                print(paste("insufficient data for prior parameter =",prior,"."),quote=FALSE)
                print("need J-q-p-1 + prior > 0", quote=FALSE)
            }
            prior<-max(prior,-(p+1))
        }
    }
    else{
        prior<--(p+1)
    }
    ## Assuming a Uniform prior on A; p(A)dA propto dA.
    ## For B0 = (I+A*)^{-1}, this implies p(B0)dB0 \propto |B0|^{-(p+1)}dB0.
    if(J-q+prior-p-1 <= 0)
        prior<-0
    if(J-q+prior-p-1 <= 0)
        prior<-p+1
    ##
    list(Y=Y,V=V,w=w,J=J,p=p,q=q,prior=prior)
}


standard.f <- function(Y, V, V0) {
        s <- svd(V0)
        n <- length(s$d)
        rtV0 <- s$u %*% diag(sqrt(s$d), n, n) %*% t(s$v)
        Ys <- solve(rtV0, Y)  ## rtV0^{-1} Y
        Vs <- array(dim = dim(V))
        for(i in seq_len(dim(Vs)[3])) {
                ## rtV0^{-1} V rtV0^{-1}
                Vs[, , i] <- solve(rtV0, t(solve(rtV0, V[, , i])))
        }
        list(Y = Ys, V = Vs, rtVo = rtV0)
}

mammult<-function(M,A,tM){
    ## Calls Fortran subroutine "mammult" to pre-multiply by M
    ## (pxp) and post-multiply by tM (pxp) each of the N pxp  
    ## elements of A (pxpxN).
    ##
    pp<-dim(A)
    p<-pp[1]
    N<-pp[3]
    MAM<-A*0.0
    out<-.Fortran("mammult",
                  as.double(M),
                  as.double(tM),
                  as.double(A),
                  as.integer(p),
                  as.integer(N),
                  as.double(0*M), ## U
                  as.double(0*M), ## V
                  MAM=MAM,
                  PACKAGE = "tlnise"
                  )
    out$MAM
}


tr<-function(A){
    ## calculate trace(A):
    sum(diag(A))
}

ldet<-function(A){
    ## calculate log(det(A)):
    sum(log(abs(diag(qr(A)[[1]]))))
}
