nparncpp=function(p,
            breaks=min(2000,round(length(p)/5)),
            test=c("t","z"),
            df,
            alternative=c("two.sided", "less", "greater"),
            compromise.n=1,
            lambdas=#if(penalty_type==1)10^seq(-2,6,length=6) else 
                    10^seq(-4,6,length=11),
            deltamax='auto',
            nknots,
            ndelta=500,
            solver=c("lsei","LowRankQP","solve.QP","ipop"),
            weights=1,
            keep.cdf=NULL,
            LowRankQP.method=c('LU','CHOL'),
            lsei.method=c('chol','svd','eigen'),
            debugging=FALSE,
            ...)
{
## these are extra-arguments not mentioned in the paper
    withpi0=TRUE
    GCVfact=1
    nskip=0
    penalty_type=2
## CHECK solver
    solver=match.arg(solver)
    if(solver=="LowRankQP") {
######## this is commented out as pi0 contains a copy of LowRankQP.c
        stopifnot(requireNamespace("LowRankQP"),quitely=TRUE)
        LowRankQP.method=match.arg(LowRankQP.method)
    }else if(solver=="solve.QP") {
        stopifnot(requireNamespace("quadprog"),quitely=TRUE)
    }else if (solver=='ipop') {
        stopifnot(requireNamespace("kernlab"),quitely=TRUE)
    }else if (solver=='lsei') {
        stopifnot(requireNamespace("limSolve"),quitely=TRUE)
        lsei.method=match.arg(lsei.method)
    }else stop("solver unimplemented!")

    alternative=match.arg(alternative);
    test=match.arg(test)
    if (test=="t" & missing(df)) stop("df is missing for t-test")
    if (test=="z") df=Inf
    if(length(breaks)==1){
            breaks=seq(0,1,length=breaks+1)
    }else{
        if(min(breaks)>0)breaks=c(0,breaks)
        if(max(breaks)<1)breaks=c(breaks,1)
    }
    counts=hist(p,breaks=breaks,plot=FALSE)$counts
    nobs=length(p)
    nbin=length(breaks)-1

    if(deltamax=='auto'){
        deltamax=0
        repeat{
            deltamax=deltamax+1
            if(cond.cdf(p.eval=breaks[2],ncp=deltamax,
                df=df,test=test,alternative=alternative,keep.cdf=NULL)>.95)
                break
        }
    }
    if(missing(nknots))nknots=max(8,2*round(deltamax))

## Begin
        binnum=1:nbin
        weights=weights*(binnum>nskip)


        binwidths=diff(breaks)
#        binwid=binwidths[1]
        bincenters=(breaks[1:nbin+1]+breaks[1:nbin])/2

        delta=seq(.Machine$double.eps,deltamax,length=ndelta)
        knots=seq(0,deltamax,length=nknots)
        tmp=NBsplines(delta,knots,1);b=tmp$b;newknows=tmp$newknots
        B=2*b[,2:(ncol(b)-3),drop=FALSE] ## b has AUC 0.5
		cdfp.tmp=cond.cdf(p.eval=breaks,ncp=delta,test=test,
                    alternative=alternative,df=df,keep.cdf=keep.cdf) #save,
        dif_cdfp=cdfp.tmp[1:nbin+1,]-cdfp.tmp[1:nbin,]
        KK=ncol(B)
        B=cbind(1,B)
        Z=matrix(1,nbin,KK+1) ## binwidths/binwidths==1
        for(j in 1:KK){
            newcol=dif_cdfp%*%B[,j+1]
            newcol=newcol/sum(newcol) /binwidths;
#            B[,j+1]=nbin*B[,j+1]/sum(newcol)	#b2(:,j+1) = nbin* b2(:,j+1) /sum(newcol) ;
            Z[,j+1]=newcol
        }

        A=diag(KK+1)
        A[2,2]=2
        ### added
        A[1,1]=0
        ### end of adding

        if(penalty_type==1){
            D=diag(0,nknots)
            for(j in 2:(nknots-1)){
                D[j,j]=1
                D[j,j+1]=-1
            }
        }else if(penalty_type==2){
            D=diag(0,nknots-3,nknots)
            for(j in 1:(nknots-3)){
                D[j,j+1]=1;
                D[j,j+2]=-2
                D[j,j+3]=1
            }
        }

        ADDA=tcrossprod(tcrossprod(A,D)) ## A%*%t(D)%*%D%*%A
        y=counts/(nobs*binwidths)

        #W=diag(weights,nbin,nbin)        ## this is not efficient in terms of space complexity; 
        #Wy=W%*%y                        ## so the original implementation is better.
        Wy=weights*y
        yWWy=crossprod(Wy)   ## Wy%*%Wy
        #WZ=W%*%Z
        WZ=weights*Z
        yWWZ=Wy%*%WZ
        ZWWZ=crossprod(WZ)   ## t(WZ)%*%WZ

        if(solver=="solve.QP"){
          Amat=cbind(rep(1,nknots),diag(1,nknots))
          bvec=rep(1:0,c(1,nknots))
          meq=ifelse(withpi0,1,2)
        }else if(withpi0 && solver=="LowRankQP"){
          Amat=matrix(1,1,nknots)
          bvec=1
          uvec=rep(1,nknots)
        }else if(!withpi0 && solver=="LowRankQP"){
          Amat=rbind(c(1,rep(0,nknots-1)),rep(1,nknots))
          bvec=c(0,1)
          uvec=rep(1,nknots)
        }else if(withpi0 && solver=="ipop"){ #withpi0
          Amat=matrix(1,1,nknots)
          bvec=1
          rvec=0
          uvec=rep(1,nknots)
          lvec=rep(0,nknots)
        }else if(!withpi0 && solver=="ipop"){ #pi0=0
          Amat=rbind(c(1,rep(0,nknots-1)),rep(1,nknots))
          bvec=c(0,1)
          rvec=c(0,0)
          uvec=rep(1,nknots)
          lvec=rep(0,nknots)
        }else if (withpi0 && solver=='lsei'){
          Emat=matrix(1,1,nknots)
          fvec=1
          Gmat=diag(nknots)
          hvec=rep(0,nknots)
        }else if (!withpi0 && solver=='lsei'){
          Emat=rbind(c(1,rep(0,nknots-1)),rep(1,nknots))
          fvec=c(0,1)
          Gmat=diag(nknots)
          hvec=rep(0,nknots)
        }

        thetahat=matrix(,length(lambdas),nknots)
        gcv=p0hat=eff.df=mins=mins.nopen=rep(NA,length(lambdas))

        for(ilam in 1:length(lambdas)){
            lambda=lambdas[ilam]

                f=-yWWZ
                H=ZWWZ+lambda*ADDA

            if(solver=="solve.QP"){
               thetanew=solve.QP(Dmat=H,
                                dvec=yWWZ,
                                Amat=Amat,
                                bvec=bvec,
                                meq=meq
                        )
               curmin=thetanew$value
               thetanew=thetanew$solution
               curmin.nopen=-yWWZ%*%thetanew+.5*t(thetanew)%*%ZWWZ%*%thetanew
            }else if(solver=="LowRankQP"){
               thetanew=as.vector(
                        LowRankQP::LowRankQP(Vmat=H,
                                  dvec=-yWWZ,
                                  Amat=Amat,
                                  bvec=bvec,
                                  uvec=uvec,
                                  method=LowRankQP.method
                        )$alpha)
               curmin=-yWWZ%*%thetanew+.5*t(thetanew)%*%H%*%thetanew
               curmin.nopen=-yWWZ%*%thetanew+.5*t(thetanew)%*%ZWWZ%*%thetanew
            }else if (solver=='ipop'){
                thetanew=primal(ipop(c=-yWWZ, 
                                     H=H, 
                                     A=Amat, 
                                     b=bvec, 
                                     l=lvec, 
                                     u=uvec, 
                                     r=rvec,
                         ...))
                curmin=-yWWZ%*%thetanew+.5*t(thetanew)%*%H%*%thetanew
                curmin.nopen=-yWWZ%*%thetanew+.5*t(thetanew)%*%ZWWZ%*%thetanew
            }else if (solver=='lsei'){
                if(lsei.method=='chol'){
                    Amat=chol(H)
                    bvec=backsolve(Amat,drop(yWWZ),transpose=TRUE)
                }else if(lsei.method=='svd'){
                    tmp=svd(H,nv=0)
                    Amat=tcrossprod(tmp$u%*%diag(sqrt(tmp$d)),tmp$u)
                    bvec=solve(Amat,drop(yWWZ))
                }else if(lsei.method=='eigen'){
                    tmp=eigen(H)
                    Amat=tcrossprod(tmp$vec%*%diag(sqrt(tmp$val)),tmp$vec)
                    bvec=solve(Amat,drop(yWWZ))
                }
                thetanew=limSolve::lsei(A=Amat, B=bvec, E=Emat, F=fvec, G=Gmat, H=hvec,...)
                curmin=thetanew$solutionNorm-crossprod(bvec)/2
                thetanew=thetanew$X
                curmin.nopen=-yWWZ%*%thetanew+.5*t(thetanew)%*%ZWWZ%*%thetanew
            }

            thetahat[ilam,]=thetanew
            mins[ilam]=curmin
            mins.nopen[ilam]=curmin.nopen
                eff.df[ilam]=sum(diag(solve(H)%*%ZWWZ))

            p0hat[ilam]=thetahat[ilam,1]

                gcv[ilam]=(yWWy-2*yWWZ%*%thetanew+thetanew%*%ZWWZ%*%thetanew)/
                                (nbin-nskip-GCVfact*eff.df[ilam])^2

        }
        imin=max(which(gcv==min(gcv)))
        attr(gcv,'which.min')=imin
        attr(gcv,'lambdas')=lambdas
        attr(gcv,'effective.df')=eff.df
        thetahatmin=thetahat[imin,]
        min=mins[imin]
        min.nopen=mins.nopen[imin]
        fhat=Z%*%thetahatmin

        compromise=Z[,1:(compromise.n+1)]%*%thetahatmin[1:(compromise.n+1)]
	compromise.min=min(compromise)
        fhat.min = min(fhat) ;
        pi0=thetahatmin[1]
        ghat=B[,2:ncol(B)]%*%thetahatmin[2:length(thetahatmin)]
        ghat=ghat/(1-thetahatmin[1])
        fhatp0 = pi0 + 0 * fhat ;
        if(1-thetahatmin[1]>0){
            mixing.prop=thetahatmin[-1]/(1-thetahatmin[1])
        }else{
            mixing.prop=thetahatmin[-1]
        }

        pdf.p=approxfun(c(0,bincenters,1),
                        c(fhat[1],fhat,tail(fhat,1)),
                        yleft=0,yright=0)
        cdf.p=approxfun(breaks,
                        c(0,cumsum(binwidths*fhat)),
                        yleft=0,yright=1)
        EDR=function(p)ifelse(p>=0 & p<=1, (cdf.p(p)-p*pi0)/(1-pi0),0)
        pdf.ncp=approxfun(c(0,delta),
                          c(ghat[1],ghat),
                          yleft=0,yright=0)
        cdf.ncp=approxfun(c(0,delta),
                          cumsum(c(0,ghat*c(0,diff(delta)))),
                          yleft=0,yright=1)
        cdf.p.p=cdf.p(p)
        pdf.p.p=pdf.p(p)
        FDR=cbind(p*pi0/cdf.p.p,
                  p*compromise.min/cdf.p.p,
                  p*fhat.min/cdf.p.p)
        LFDR=cbind(pi0/pdf.p.p,compromise.min/pdf.p.p,fhat.min/pdf.p.p)
        TN=cbind(pi0*(1-p)/(1-cdf.p.p),
                 compromise.min*(1-p)/(1-cdf.p.p),
                 fhat.min*(1-p)/(1-cdf.p.p))
        colnames(FDR)=colnames(LFDR)=colnames(TN)=c('pi0','compromise','f1')

        FIDR=function(delta.interest){
            idx=delta<=delta.interest
            numerator=breaks*pi0+(1-pi0)*
               if(sum(idx)>0)(cdfp.tmp[,idx]%*%ghat[idx])*(deltamax/ndelta)else 0
            numerator=approx(breaks,c(numerator),p)$y
            pmax(0,pmin(1,numerator/cdf.p.p))
        }


        rslt=list(pi0=pi0,
                  compromise=compromise.min,
                  f1=fhat.min,
                  pdf.p=pdf.p,
                  cdf.p=cdf.p,
                  EDR=EDR,
                  pdf.ncp=pdf.ncp,
                  cdf.ncp=cdf.ncp,
                  FDR=FDR,
                  LFDR=LFDR,
                  FIDR=FIDR,
                  TP=1-FDR,
                  TN=TN,
                  mixing.prop=mixing.prop,
                  agcv=gcv,
                  df=df,
                  test=test,
                  alternative=alternative,
                lambdas=lambdas,
                ndelta=ndelta,
                deltamax=deltamax,
                nknots=nknots,
                bincenters=bincenters,
                weights=weights,
                solver=solver, 
                par=thetahat, 
                data=list(p=p),
                nobs=nobs
             )
        if(isTRUE(debugging))rslt=c(rslt,list(
                  thetahat=thetahat,
                  thetahatmin=thetahatmin,
                  fhat=fhat,
                  fhat.min=fhat.min,
                  p0hat=p0hat,
                  imin=imin,
                  ZWWZ=ZWWZ,
                  yWWZ=yWWZ,
                  KK=KK,
                  D=D,
                pi0=pi0,
                b=b,
                knots=knots,
                counts=counts,
                binwidths=binwidths,
                counts=counts,
                ghat=ghat,
                B=B,
                GCVfact=GCVfact,
                ADDA=ADDA,
                A=A
                ,Z=Z
                ,compromise.min=compromise.min,
                compromise=compromise,
                fhatp0=fhatp0,
                cdfp=cdfp.tmp,
                b=b,
                withpi0=withpi0,
                min=min,
                min.nopen=min.nopen
        ))
        class(rslt)=c('nparncpp',"ncpest")
        rslt
}
