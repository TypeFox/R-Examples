sparncpF=function(obj1, obj2, ...) UseMethod('sparncpF')

sparncpF.nparncpF=function(obj1, obj2,...)
{   tmp=obj1;
    obj1=obj2;
    obj2=tmp
    sparncpF.parncpF(obj1, obj2, ...)
}

sparncpF.numeric=function(obj1, obj2, ...)
{   Fstat=obj1; df1=obj2
    obj1=parncpF(Fstat,df1, central=FALSE, ...)
    obj2=nparncpF(Fstat,df1, ...)
    sparncpF.parncpF(obj1, obj2, ...)
}

sparncpF.parncpF=function(obj1, obj2, ...)
{   parfit=obj1; nparfit=obj2
    if(!inherits(nparfit,'nparncpF')) stop("second argument should be an 'nparncpF' object")
    if(!all.equal(parfit$data, nparfit$data)) stop("the two objects are based on different data")

#    obj=function(propar) propar*marginal.df(parfit)(parfit$data$Fstat)+(1-propar)*marginal.df(nparfit)(nparfit$data$Fstat) ## this uses mixture density for different df's
    fitted.parfit=fitted(parfit); fitted.nparfit=fitted(nparfit)
    obj=function(propar) sum( log(propar*pmax(fitted.parfit,0)+(1-propar)*pmax(0,fitted.nparfit) ))
    propar.fit=optimize(obj, c(0,1),  maximum =TRUE,tol=5e-4)
    propar=round(propar.fit$maximum,3)
    ll=propar.fit$objective
    attr(ll,'df')=ifelse(propar==1 || propar==0, 0,1)+propar*parfit$enp+(1-propar)*nparfit$enp
	attr(ll,'nobs')=nobs(obj1)
	class(ll)='logLik'

    ans=list(pi0=propar*parfit$pi0+(1-propar)*nparfit$pi0,
             mu.ncp=propar*parfit$mu.ncp+(1-propar)*nparfit$mu.ncp,
             sd.ncp=sqrt((propar*parfit$sd.ncp)^2+((1-propar)*nparfit$sd.ncp)^2),
             logLik=ll, enp=attr(ll,'df'), par=propar, 
             gradient=NULL,      hessian=NULL,      ## TO BE ADDED LATER
             parfit=parfit, nparfit=nparfit, nobs=nobs(obj1)
             )
    class(ans)=c('sparncpF','ncpest')
    ans
}

fitted.sparncpF=#fitted.values.sparncpF=
function(object, ...)
{
    object$par * fitted(object$parfit) + (1-object$par)*fitted(object$nparfit)
}
print.sparncpF=function(x,...)
{
    cat('pi0=',x$pi0, fill=TRUE)
    cat('mu.ncp=', x$mu.ncp, fill=TRUE)
    cat('sd.ncp=', x$sd.ncp, fill=TRUE)
    cat("rho=", x$par, fill=TRUE)
    cat("enp=", x$enp, fill=TRUE)
#    cat("lambda=", x$lambda, fill=TRUE)
    invisible(x)
}
summary.sparncpF=function(object,...)
{
    print.sparncpF(object,...)
}
plot.sparncpF=function(x,...)
{
#    dev.new(width=8, height=4)
    op=par(mfrow=c(1,2))
    hist(x$parfit$data$Fstat, pr=TRUE, br=min(c(max(c(20, length(x$parfit$data$Fstat)/100)), 200)), xlab='t',main='t-statistics')
    ord=order(x$parfit$data$Fstat)
    lines(x$parfit$data$Fstat[ord], fitted(x)[ord], col='red', lwd=2)
    d.ncp=function(d) dncp(x)(d)
    curve(d.ncp, quantile(x$parfit$data$Fstat,.001), quantile(x$parfit$data$Fstat,.999), 500, 
        xlab=expression(delta), ylab='density',main='noncentrality parameters')
    abline(v=c(0, x$mu.ncp), lty=1:2)
    par(op)
    invisible(x)
}

#sparncpF.numeric=function(Fstat, df, parfit, ...)
#{
#    method='SQP'
#    if (method=='SQP') {
#        sparncpF.sqp(Fstat, df, parfit, ...)
#    }
#}

#sparncpF.sqp = function (Fstat, df, parfit=parncpF(Fstat,df,zeromean=F), penalty=c('3rd.deriv','2nd.deriv','1st.deriv'), lambdas=10^seq(-1,15,length=17), starts, smooth.enp=FALSE, IC=c('BIC','CAIC','HQIC','AIC'),
#                        K=100, bounds=quantile(Fstat,c(.01,.99)), solver=c('solve.QP','lsei','ipop','LowRankQP'),plotit=FALSE, verbose=FALSE, approx.hess=TRUE, ... )
#{
##   source("int.nct.R"); source("laplace.nct.R"); source("saddlepoint.nct.R"); source("dtn.mix.R")
#    solver=match.arg(solver)
#    penalty=match.arg(penalty)
#    solver.package=switch(solver, solve.QP='limSolve', ipop='kernlab', lsei='limSolve',LowRankQP='LowRankQP')
#library(solver.package)
#library("Matrix")
#    if (K<=0 || length(K)!=1) stop("K should be a positive integer")
#    if (any(lambdas<0)) stop("lambdas should be a vector of positive numbers")
#    if (length(bounds)==1) bounds=c(-abs(bounds), abs(bounds))
#    if (diff(bounds)[1]<=0) bounds=bounds[2:1]
#
#    G=max(c(length(Fstat),length(df)))
#    mus=seq(bounds[1], bounds[2], length=K)
#    sigs=rep(diff(bounds)[1]/(K-1), K)
#
#    IC=toupper(IC); IC=match.arg(IC)
#    IC.fact=switch(IC, AIC=2, BIC=log(G), CAIC=log(G)+2, HQIC=2*log(log(G)))
#
#    h0.i=df(Fstat,df) 
#    b.i=matrix(0,G,K)
#    for(k in 1:K){
#        b.i[,k]=dtn.mix(Fstat,df,mus[k], sigs[k],...)
#        if(any(b.i[,k]<0) || any(is.na(b.i[,k]))) {
#            warning("Noncentral density unreliable. I switched to exact density function")
#            b.i[,k]=dtn.mix(Fstat,df,mus[k], sigs[k],approximation='none')
#        }
#    }
#    if(approx.hess==1) {
#        b.approx=as(b.i*(b.i>mean(b.i)), 'dgCMatrix')
#    }else if(approx.hess>0 && approx.hess<1) {
#        b.approx=as(b.i*(b.i>quantile(b.i,approx.hess)),'dgCMatrix')
#    }
#    b.i=b.i-h0.i
#
#    if(penalty=='1st.deriv'){
##        Omega=diag(0,K,K)                   ### 1st order derivative
##        for(j in 1:K) for(k in 1:j) {
##                tmp=(mus[j]-mus[k])^2/(sigs[j]^2+sigs[k]^2); 
##                Omega[j,k]=Omega[k,j]=(1-tmp)/sqrt(2*pi*(sigs[j]^2+sigs[k]^2)^3)*exp(-tmp/2)
##        }
#            diffmu=outer(mus,mus,'-'); sumsig2=outer(sigs^2, sigs^2, '+')
#            Omega=exp(-diffmu^2/2/sumsig2)/sqrt(2*pi*sumsig2^5)*(sumsig2-diffmu^2)
#        tmp=(mus-parfit$mu)^2/(sigs^2+parfit$sd^2)
#        omega=(1-tmp)/sqrt(2*pi*(sigs^2+parfit$sd^2)^3)*exp(-tmp/2)
#        constant=1/4/sqrt(pi)/parfit$sd^3
#    }else if(penalty=='2nd.deriv'){
#            diffmu=outer(mus,mus,'-'); sumsig2=outer(sigs^2, sigs^2, '+')
#            Omega=exp(-diffmu^2/2/sumsig2)/sqrt(2*pi*sumsig2^9)*(3*sumsig2^2-6*sumsig2*diffmu^2+diffmu^4)
#        diffmu2=(mus-parfit$mu)^2; sumsig2=(sigs^2+parfit$sd^2)
#        omega=exp(-diffmu2/2/sumsig2)/sqrt(2*pi*sumsig2^9)*(3*sumsig2^2-6*sumsig2*diffmu2+diffmu2^2)
#        constant=3/8/sqrt(pi)/parfit$sd^5
#    }else if(penalty=='3rd.deriv'){
#            diffmu=outer(mus,mus,'-'); sumsig2=outer(sigs^2, sigs^2, '+')
#            Omega=exp(-diffmu^2/2/sumsig2)/sqrt(2*pi*sumsig2^13)*(15*sumsig2^3-45*sumsig2^2*diffmu^2+15*sumsig2*diffmu^4-diffmu^6)
#        diffmu2=(mus-parfit$mu)^2; sumsig2=(sigs^2+parfit$sd^2)
#        omega=exp(-diffmu2/2/sumsig2)/sqrt(2*pi*sumsig2^13)*(15*sumsig2^3-45*sumsig2^2*diffmu2+15*sumsig2*diffmu2^2-6*diffmu2^3)
#        constant=15/16/sqrt(pi)/parfit$sd^7
#    }
#
#    NPLL=function(thetas, take.sum=TRUE ){ ### depends on h0.i, b.i, lambda, Omega, G, omega
#        one_pi0=max(c(sum(thetas), 1e-6))
#        Lik=pmax(h0.i+b.i%*%thetas, 1e-8)  #### CHECK negativity PROBLEM
#        NLL=-log(Lik)
#        if(take.sum) ans=sum(NLL)+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)-lambda/one_pi0*drop(crossprod(omega,thetas))+.5*lambda*constant
#        else ans=NLL+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)/G-lambda/one_pi0*drop(crossprod(omega,thetas))/G+.5*lambda*constant/G
#        ans
#    }
#
#    grad.NPLL=function(thetas, take.sum=TRUE) {    ### depends on h0.i, b.i, lambda, Omega, G, K
#        one_pi0=max(c(sum(thetas), 1e-6))
#        Lik=drop(h0.i+b.i%*%thetas)     #### no log is taken later
#        Omega.theta=drop(Omega%*%thetas)
#        if(take.sum){
#           -drop(crossprod(1/Lik, b.i))+
#           lambda/one_pi0/one_pi0*Omega.theta - 
#           lambda*drop(crossprod(thetas,Omega.theta))/one_pi0/one_pi0/one_pi0 -
#           lambda/one_pi0*omega + lambda*drop(crossprod(omega, thetas))/one_pi0/one_pi0
#        }else{  # return a G x k matrix, each row is the gradient based on one data point
#            -b.i/Lik+ 
#            matrix(lambda/one_pi0/one_pi0*Omega.theta/G, G ,K, byrow=TRUE)-
#            lambda*drop(crossprod(thetas,Omega.theta))/one_pi0/one_pi0/one_pi0/G -
#            matrix(lambda/one_pi0*omega/G, G, K, byrow=TRUE) +
#            lambda*drop(crossprod(omega,thetas))/one_pi0/one_pi0/G
#        }    
#    }
#        
#    hess.NPLL=function(thetas,approx=FALSE){         ### depends on h0.i, b.i, lambda, Omega, K
#        one_pi0=max(c(sum(thetas),1e-6))        ## avoiding division by zero
#        Lik= drop(h0.i+b.i%*%thetas)
#        Omega.theta=drop(Omega%*%thetas)
#        Omega.theta.onep=matrix(Omega.theta,K,K)
#        omega.onep=matrix(omega, K, K)
#        if(approx){
#            bl.approx=b.approx/Lik
#            h0l=h0.i/Lik
#            blh0l=matrix(as.vector(crossprod(bl.approx, h0l)),K,K)
#            crossprod.term=as.matrix(crossprod(bl.approx))+drop(crossprod(h0l))-blh0l-t(blh0l)
#        }else{
#            bi.lik=b.i/Lik
#            crossprod.term= crossprod(bi.lik)
#        }
#        crossprod.term + lambda*Omega/one_pi0/one_pi0 - 
#            2*lambda*(Omega.theta.onep+t(Omega.theta.onep))/one_pi0/one_pi0/one_pi0+
#            3*lambda*matrix(drop(crossprod(thetas,Omega.theta)),K,K)/one_pi0/one_pi0/one_pi0/one_pi0 +
#            lambda/one_pi0/one_pi0*(omega.onep+t(omega.onep)) +
#            2*lambda*drop(crossprod(omega, thetas))/one_pi0/one_pi0/one_pi0
#    }
#
#    C.fctn=function(thetas, eps=1e-8){ ## constraint function s.t. C.fctn(thetas)>=0
#        rbind(diag(1,length(thetas)),-1)%*%thetas+rep(c(eps,1-eps),c(length(thetas),1))
#    }
#    grad.C=function(thetas){  ## grad.C^TRUE thetas + C >=0
#        cbind(diag(1,length(thetas)),-1)
#    }
#    Amat=grad.C(numeric(K)) ## this is the A matrix for quadprog::solve.QP, i.e., t(A)%*%theta>=theta0 linear constraints
#                            ## for limSolve::lsei(, this is t(G)
#
#
#    sqp=function(thetas, conv.f=1e-10, verbose=FALSE, maxiter=1e3) { ## thetas is a starting value
#        ## not a general solver; instead, designed for this problem per se
#        ## depends on solver, grad.NPLL, hess.NPLL, nearPD, Amat, C.fctn
#
#        npll.last=Inf
#        niter=1
#        repeat{
#            thetas=pmax(thetas,0)
#            dvec=-grad.NPLL(thetas)    ## negative gradiant
#            fnscale=10^floor(max(log10(abs(dvec))))
#            dvec=dvec/fnscale
#        #    tmp=drop(1/((1-sum(betas))*F0+F1%*%betas))*(F1-F0)
#            Dmat=hess.NPLL(thetas, approx=approx.hess>0)/fnscale
#        #    repeat{
#        #        eigv=eigen(Dmat,TRUE)
#        #        if(tail(eigv$val,1)<1e-5) {
#        #            if(verbose)cat('\tmin.eigval=',tail(eigv$val,1),'\trk.D=',qr(Dmat)$rank,fill=TRUE)
#        #            Dmat=eigv$vec%*%diag(eigv$val-min(eigv$val)+1.1e-6)%*%t(eigv$vec) 
#        #        }else 
#        #            break
#        #    }
#            Dmat=as.matrix(nearPD((Dmat+t(Dmat))/2 )$mat)
#        #    Amat=t(UI)
#        #    bvec=drop(-UI%*%thetas+CI) -1e-10
#            bvec=-C.fctn(thetas)
#            qp.sol=if(solver=='solve.QP') {
#                        tmpA=try(solve.QP(Dmat,dvec,Amat=Amat,bvec)$solution, silent=TRUE)
#                        if(class(tmpA)=='try-error'){ # solver='lsei'
#                            tmpA=chol(Dmat); tmpB=.5*drop(forwardsolve(t(tmpA), dvec)); 
#                            limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,0,K), F=numeric(0), G=t(Amat), H=bvec, 
#                                 Wx=NULL, Wa=NULL, type=1)$X
#                        }else tmpA
#                   }else if (solver=='lsei') {
#                        tmpA=chol(Dmat); tmpB=.5*drop(forwardsolve(t(tmpA), dvec)); 
##                        limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,1,K), F=0, G=t(Amat), H=bvec, 
##                             Wx=NULL, Wa=NULL, type=1)$X
#                        limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,0,K), F=numeric(0), G=t(Amat), H=bvec, 
#                             Wx=NULL, Wa=NULL, type=1)$X
#                   }else if (solver=='ipop') {  ## not working well ## the R translation of the LOQO code is not very honest
#                        tmpA=ipop(c=-dvec, H=Dmat, A=t(Amat), b=bvec, l=rep(-1,K), u=rep(1,K), r=rep(1e6,K+1),verb=verbose)
#                        if(tmpA@how!='converged'){warning('ipop not converged')}
#                        tmpA@primal
#                   }else if (solver=='LowRankQP') { ## not working well
#                        tmpTransform=solve(tcrossprod(Amat),Amat)
#                        Hmat=crossprod(tmpTransform, Dmat%*%tmpTransform)
#                        tmpAns=LowRankQP(Vmat=Hmat,dvec=Hmat%*%bvec-crossprod(tmpTransform,dvec),
#                                        Amat=matrix(0,0,K+1),bvec=numeric(0),uvec=rep(1e6,K+1),method='LU')
#                        drop(tmpTransform%*%(bvec+tmpAns$alpha))
#                   }
#            thetas.new=thetas+qp.sol
#            repeat{
#                if(all(C.fctn(thetas.new)>=0-1e-10) )break  #&& crossprod(Amat,thetas.new)>=bvec
#                thetas.new=(thetas+thetas.new)/2
#                if(verbose) cat("halving due to constraints",fill=TRUE)
#            }
#            nhalving=1
#            repeat{
#                npll.new=NPLL(thetas.new)
#                if(npll.new<=npll.last) break
#                if(all(thetas==thetas.new)) break
#                thetas.new=(thetas+thetas.new)/2
#                if(verbose) cat("halving due to increase in objective",fill=TRUE)
#                nhalving=nhalving+1
#                if(nhalving>50) {warning("step-halving limite reached"); break}
#            }
#            thetas=thetas.new
#            if(verbose)cat('npll.new=', npll.new, '\tpi0=',1-sum(thetas), '\tmin.theta=', min(thetas), '\tfnscale=', fnscale, fill=TRUE)
#            if(npll.last-npll.new<conv.f*fnscale) break
#            niter=niter+1
#            if(niter>maxiter){warning("maximum iteration reached; results from last iteration returned"); break}
#            npll.last=npll.new
#         }
#         if(verbose) cat("DONE!",fill=TRUE)
#         return(pmax(thetas,0))
#    }
#
#    enp=function(thetas, eps=1e-6)  ## effective number of parameters, ## depends on hess.NPLL, grad.NPLL
#    {
#        nonzero.idx=which(thetas>eps)
#
#        if(length(nonzero.idx)==0){        # pi0=1
#            return(0)
#        }
#        hess=hess.NPLL(thetas)[nonzero.idx,nonzero.idx,drop=FALSE]
#        grads=grad.NPLL(thetas,take.sum=FALSE)[,nonzero.idx,drop=FALSE]
#        Kmat=crossprod(grads)
#
#        if( 1-sum(thetas)<eps) { ## pi0=0, i.e. sum(thetas)==1; ## this needs reparameterization
#            if(length(nonzero.idx)==1) return(0) ## only one beta is nonzero, but pi0=0 too
#            Jacobian=rbind(-1,diag(1,length(thetas[nonzero.idx])-1))
#            hess=t(Jacobian)%*%hess%*%Jacobian
#            Kmat=t(Jacobian)%*%Kmat%*%Jacobian
#        }
#
#        return(max(c(0, sum(diag(solve(nearPD(hess)$mat, Kmat)))+2-(parfit$mu==0))))
#    }
#
#    if(missing(starts)) {
#        tmpobj=function(fact) abs(fact*parfit$mu-sum(mus*dnorm(mus,parfit$mu,parfit$sd)))
#        fact=optimize(tmpobj,c(1e-3,1e2))$minimum
#        beta=dnorm(mus,parfit$mu,parfit$sd)/fact
#        beta=beta/sum(beta)
#        starts=(1-min(c(.99,max(c(.01,parfit$pi0)))))*beta
#    }
#
#    n.lambda=length(lambdas)
#    nics=enps=nic.sd=pi0s=numeric(n.lambda)
#    sqp.fit=matrix(NA_real_, n.lambda, K)
#    for(i in n.lambda:1) {
#        lambda=lambdas[i]
#        sqp.fit[i,]=sqp(if(i==n.lambda) starts else sqp.fit[i+1,], conv.f=1e-6, verbose=verbose)
#        enps[i]=enp(sqp.fit[i,])
#        tmp=NPLL(sqp.fit[i,],FALSE)
#        nics[i]=2*(sum(tmp))+IC.fact*enps[i]
#        nic.sd[i]=sd(tmp)*sqrt(G)
#        pi0s[i]=1-sum(sqp.fit[i,])
#    }
#
#    if(smooth.enp) {
#          warning("smoothing snp is not well tested")
#library("monoProc")
#        loe=loess(enps~log10(lambdas))
#        mon=mono.1d(list(log10(lambdas), fitted(loe)), bw.nrd0(fitted(loe))/3,mono1='decreasing')
#        enps.smooth=mon@y
#        nics=nics-enps+enps.smooth
#    }
#
#    i.final=which.min(nics); # i.1se=tail(which(nics<=nics[i.final]+nic.sd[i.final] & 1:n.lambda>=i.final),1) # one se rule
#    lambda.final=lambdas[i.final]
#    beta.final=sqp.fit[i.final,]/sum(sqp.fit[i.final,])
#
#
#    ll=-NPLL(sqp.fit[i.final,])
#    attr(ll, 'df')=enps[i.final]
#    class(ll)='logLik'
#
#    ans=list(pi0=pi0s[i.final], mu.ncp=beta.final%*%mus, sd.ncp= sqrt(beta.final%*%(mus^2+sigs^2)-(beta.final%*%mus)^2), 
#             logLik=ll, enp=enps[i.final], par=sqp.fit[i.final,],lambda=lambdas[i.final],
#             gradiant=grad.NPLL(sqp.fit[i.final,]), hessian=hess.NPLL(sqp.fit[i.final,]),
#
#             beta=beta.final, IC=IC, 
#             all.mus=mus, all.sigs=sigs, data=list(Fstat=Fstat, df=df), i.final=i.final, all.pi0s=pi0s,
#             all.enps=enps, all.thetas=sqp.fit, all.nics=nics, all.nic.sd=nic.sd, all.lambdas=lambdas)
#    class(ans)=c('sparncpF','ncpest')
#    if(plotit) plot.sparncpF(ans)
#    ans
#}

#vcov.ncpest=function(object,...)
#{
#    object$hessian
#}
#logLik.ncpest=function(object,...)
#{
#    object$logLik
#}
#coef.ncpest=#coefficients.ncpest=
#function(object,...)
#{
#    object$par
#}
#fitted.sparncpF=#fitted.values.sparncpF=
#function(object, ...)
#{
#    
#    nonnull.mat=matrix(NA_real_, length(object$data$Fstat), length(object$all.mus))
#    for(k in 1:length(object$all.mus)) nonnull.mat[,k]=dtn.mix(object$data$Fstat, object$data$df, object$all.mus[k],object$all.sigs[k],FALSE,...)
#    object$pi0*df(object$data$Fstat, object$data$df)+(1-object$pi0)*drop(nonnull.mat%*%object$beta)
#}

#summary.sparncpF=function(object,...)
#{
#    print.sparncpF(object,...)
#}
#print.sparncpF=function(x,...)
#{
#    cat('pi0=',x$pi0, fill=TRUE)
#    cat('mu.ncp=', x$mu.ncp, fill=TRUE)
#    cat('sd.ncp=', x$sd.ncp, fill=TRUE)
#    cat("enp=", x$enp, fill=TRUE)
#    cat("lambda=", x$lambda, fill=TRUE)
#    invisible(x)
#}
#plot.sparncpF=function(x,...)
#{
#    dev.new(width=7,heigh=7)
#    op=par(mfrow=c(2,2))
##    attach(x)
#    n.lambda=length(x$all.lambdas)
#    i.1se=tail(which(x$all.nics<=x$all.nics[x$i.final]+x$all.nic.sd[x$i.final] & 1:n.lambda>=x$i.final),1) ## one se rule
#    plot(log10(x$all.lambdas), (x$all.nics), 
#        ylim=range(c(x$all.nics+x$all.nic.sd, x$all.nics-x$all.nic.sd)[c(1:min(i.1se,x$i.final+1,n.lambda),n.lambda+1:min(i.1se,x$i.final+1,n.lambda))]), 
#        xlab='log10(lambda)', ylab='NIC')
#    for(i in 1:n.lambda) lines(rep(log10(x$all.lambdas)[i],2), x$all.nics[i]+c(-1,1)*x$all.nic.sd[i], lwd=3)
#    abline(v=log10(x$all.lambdas[x$i.final]), col=2);
#    abline(h=x$all.nics[x$i.final]+c(-1,1)*x$all.nic.sd[x$i.final], col=4)
##    abline(v=log10(x$all.lambdas[i.1se]), col=2,lty=2)
#
##    plot(log10(x$all.lambdas), x$all.enps); abline(v=log10(x$all.lambdas[c(x$i.final,i.1se)]), xlab='log10(lambda)', ylab='ENP', lty=1:2)
##    plot(log10(x$all.lambdas), x$all.pi0s); abline(v=log10(x$all.lambdas[c(x$i.final,i.1se)]), xlab='log10(lambda)', ylab='pi0', lty=1:2)
#    plot(log10(x$all.lambdas), x$all.enps,xlab='log10(lambda)',ylab='effective # parameters'); 
#        abline(v=log10(x$all.lambdas[c(x$i.final)]), xlab='log10(lambda)', ylab='ENP', lty=1)
#    plot(log10(x$all.lambdas), x$all.pi0s, xlab='log10(lambda)', ylab='pi0'); 
#        abline(v=log10(x$all.lambdas[c(x$i.final)]), xlab='log10(lambda)', ylab='pi0', lty=1)
#
#    d.ncp=function(xx)        # p in the paper
#    {   ## depends on mus, sigs
#        xx=outer(xx, x$all.mus, '-')
#        xx=sweep(xx, 2, x$all.sigs, '/')
#        d=sweep(dnorm(xx),2,x$all.sigs,'/')
#        drop(d%*%x$beta)
#    }
#    curve(d.ncp, min(x$data$Fstat),max(x$data$Fstat),100,col=4,lwd=2, xlab='delta',ylab='density')
##    detach(x)
#    rug(x$all.mus)
#    par(op)
#    invisible(x)
#}

if(FALSE) {
(pfit=parncpF(Fstat,df,FALSE))
npfit.mean=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,approx.hess=TRUE)
npfit.0=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,approx.hess=0)
npfit.75=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,plotit=FALSE,approx.hess=.75)
npfit.80=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,plotit=FALSE,approx.hess=.80)
npfit.85=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,plotit=FALSE,approx.hess=.85)
npfit.90=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,plotit=FALSE,approx.hess=.90)
npfit.95=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),verbose=FALSE,plotit=FALSE,approx.hess=.95)

npfit.mean=sparncpF(Fstat,df,starts=rep((1-pfit['pi0'])/K,K),penalty='2',verbose=FALSE,approx.hess=TRUE)

############ below is old code
    df=8
    G=20000
    pi0=.6
    G0=round(G*pi0); G1=G-G0
#    ncps=c(rbeta(floor(G1/2),2,1)*3-3, rbeta(ceiling(G1/2),1,2)*3)+5
#    ncps=rnorm(G1)
    ncps=rnorm(G1,3)    
#    ncps=rnorm(G1,c(2,-2))
    Fstat=c(rt(G0,df),rt(G1,df,ncps))

Fstat.ord=order(Fstat)
K=100
mus=c(seq(min(Fstat),max(Fstat),length=K))
#    mus=c(seq(-23,23,length=K))
sigs=c(rep(mus[3]-mus[2], K))

source("int.nct.R"); source("laplace.nct.R"); source("saddlepoint.nct.R"); source("dtn.mix.R")
h0.i=df(Fstat,df) #quantile(Fstat,.05)),
b.i=matrix(0,G,K)
time0=proc.time()
for(k in 1:K) {
#           b.i[,i]=d(Fstat,df,mus[k], sigs[k], approximation='saddlepoint', normalize='int') ## old
######### new
#    b.i[,k]=df(Fstat/sigs[k],df,mus[k]/sigs[k])/sigs[k]
#    b.i[,k]=dt.sad(Fstat/sigs[k],df,mus[k]/sigs[k], normalize='approx')/sigs[k]
    b.i[,k]=dt.lap(Fstat/sigs[k],df,mus[k]/sigs[k])/sigs[k]
}
b.i=b.i-h0.i
proc.time()-time0


#Omega=diag(0,K,K)                   ### no penalty
#D=diag(1,K,K); diag(D[,-1])=-1;     ### 1st order differencing
#    Omega=crossprod(D)
#Omega=diag(6,K,K)                   ### 2nd order differencing
#    diag(Omega[-1,])=diag(Omega[,-1])=-4
#    diag(Omega[-(1:2),])=diag(Omega[,-(1:2)])=1
Omega=diag(0,K,K)                   ### 1st order derivative
for(j in 1:K) for(k in 1:j) {
        tmp=(mus[j]-mus[k])^2/(sigs[j]^2+sigs[k]^2); 
        Omega[j,k]=Omega[k,j]=(1-tmp)/sqrt(2*pi*(sigs[j]^2+sigs[k]^2)^3)*exp(-tmp/2)
}


d.ncp=function(x, betas)        # p in the paper
{
    x=outer(x, mus, '-')
    x=sweep(x, 2, sigs, '/')
    d=sweep(dnorm(x),2,sigs,'/')
    drop(d%*%betas)
}
d.t=function(thetas)  #, pi0.ub=1-1e-6  # f in the paper
{
    pi0=1-sum(thetas); #if(pi0<pi0.lb) {pi0=pi0.lb; thetas=thetas*(1-pi0.lb)}
    drop(h0.i+b.i%*%thetas)
}

NPLL=function(thetas, take.sum=TRUE ){ #,pi0.ub=1-1e-6  #if(min(betas)<0 || max(betas)>1) return(1e9) ;ps=betas/sum(betas)
    one_pi0=sum(thetas); #if(pi0<pi0.lb) {pi0=pi0.lb; thetas=thetas*(1-pi0.lb)}
    Lik=pmax(h0.i+b.i%*%thetas, 1e-8)
    NLL=-log(Lik)
    if(take.sum) ans=sum(NLL)+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)
    else ans=NLL+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)/G
    ans
}

grad.NPLL=function(thetas, take.sum=TRUE) {    #, pi0.ub=1-1e-6
    one_pi0=sum(thetas); #if(pi0<pi0.lb) {pi0=pi0.lb; thetas=thetas*(1-pi0.lb)}
    Lik=drop(h0.i+b.i%*%thetas)
    Omega.theta=drop(Omega%*%thetas)
    if(take.sum){
       -drop(crossprod(1/Lik, b.i))+
       lambda/one_pi0/one_pi0*Omega.theta - 
       lambda*drop(crossprod(thetas,Omega.theta))/one_pi0/one_pi0/one_pi0
    }else{  # return a G x k matrix, each row is the gradient based on one data point
        -b.i/Lik+ 
        matrix(lambda/one_pi0/one_pi0*Omega.theta/G, G ,K, byrow=TRUE)-
        lambda*drop(crossprod(thetas,Omega.theta))/one_pi0/one_pi0/one_pi0/G
    }    
}
    
hess.NPLL=function(thetas){
    one_pi0=sum(thetas)
    Lik= drop(h0.i+b.i%*%thetas)
    Omega.theta=drop(Omega%*%thetas)
    Omega.theta.onep=matrix(Omega.theta,K,K)
    bi.lik=b.i/Lik
   
    crossprod(bi.lik)+lambda*Omega/one_pi0/one_pi0 - 
        2*lambda*(Omega.theta.onep+t(Omega.theta.onep))/one_pi0/one_pi0/one_pi0+
        3*lambda*matrix(drop(crossprod(thetas,Omega.theta)),K,K)/one_pi0/one_pi0/one_pi0/one_pi0
}

C.fctn=function(thetas, eps=1e-8){ ## constraint function s.t. C.fctn(thetas)>=0
    rbind(diag(1,length(thetas)),-1)%*%thetas+rep(c(eps,1-eps),c(length(thetas),1))
}
grad.C=function(thetas){  ## grad.C^TRUE thetas + C >=0
    cbind(diag(1,length(thetas)),-1)
}
Amat=grad.C(numeric(K))
    
sqp=function(thetas, conv.f=1e-10, fnscale, verbose=TRUE) {
#library("quadprog")
  npll.last=Inf
  repeat{

    dvec=-grad.NPLL(thetas)    ## negative gradiant
      if(missing(fnscale)) fnscale=10^floor(max(log10(abs(dvec))))
    dvec=dvec/fnscale
#    tmp=drop(1/((1-sum(betas))*F0+F1%*%betas))*(F1-F0)
    Dmat=hess.NPLL(thetas)/fnscale
#    repeat{
#        eigv=eigen(Dmat,TRUE)
#        if(tail(eigv$val,1)<1e-5) {
#            if(verbose)cat('\tmin.eigval=',tail(eigv$val,1),'\trk.D=',qr(Dmat)$rank,fill=TRUE)
#            Dmat=eigv$vec%*%diag(eigv$val-min(eigv$val)+1.1e-6)%*%t(eigv$vec) 
#        }else 
#            break
#    }
    Dmat=as.matrix(nearPD(Dmat )$mat)
#    Amat=t(UI)
#    bvec=drop(-UI%*%thetas+CI) -1e-10
    bvec=-C.fctn(thetas)
    qp.sol=solve.QP(Dmat,dvec,Amat=Amat,bvec)
    thetas.new=thetas+qp.sol$solution
    repeat{
        if(all(C.fctn(thetas.new)>=0-1e-10) )break  #&& crossprod(Amat,thetas.new)>=bvec
        thetas.new=(thetas+thetas.new)/2
        cat("halving due to constraints",fill=TRUE)
    }
    repeat{
        npll.new=NPLL(thetas.new)
        if(npll.new<=npll.last) break
        thetas.new=(thetas+thetas.new)/2
        cat("halving due to increase in objective",fill=TRUE)
    }
    thetas=thetas.new
    if(verbose)cat('npll.new=', npll.new, '\tpi0=',1-sum(thetas), '\tmin.theta=', min(thetas), '\tfnscale=', fnscale, fill=TRUE)
    if(npll.last-npll.new<conv.f*fnscale) break
    npll.last=npll.new
 }
 if(verbose) cat("DONE!",fill=TRUE)
 return(thetas)
}
enp=function(thetas, eps=1e-8)
{
    nonzero.idx=which(thetas>eps)
#    mu.star=mus[nonzero.idx]
#    sig.start=sigs[nonzero.idx]
#
    if(length(nonzero.idx)==0){        # pi0=1
        return(0)
    }
#    else if(1-sum(thetas)<eps) {   # pi0=0, i.e. sum(thetas)==1
#        F0.star=b.i[,nonzero.idx][1]+h0.i; 
#        nonzero.idx[which(nonzero.idx)[1]]=FALSE 
#    }else {
#        F0.star=h0.i
#    }
#    thetas.star=thetas[nonzero.idx]
#    F1.star=b.i[,nonzero.idx]+h0.i
#    Omega.star=Omega[nonzero.idx,nonzero.idx]
#
#    FW=drop(1/((1-sum(thetas.star))*F0.star+F1.star%*%thetas.star))*(F1.star-F0.star)  ## G x K
#    FWWF=crossprod(FW)
#    hess=FWWF+lambda*Omega.star
#    K=FWWF-crossprod(lambda*thetas.star*Omega.star)/G
    hess=hess.NPLL(thetas)[nonzero.idx,nonzero.idx]
    grads=grad.NPLL(thetas,take.sum=FALSE)[,nonzero.idx]
    K=crossprod(grads)

if( 1-sum(thetas)<eps) {## pi0=0, i.e. sum(thetas)==1; ## this needs reparameterization
    if(length(nonzero.idx)==1) return(0) ## only one beta is nonzero, but pi0=0 too
    Jacobian=rbind(-1,diag(1,length(thetas[nonzero.idx])-1))
    hess=t(Jacobian)%*%hess%*%Jacobian
    K=t(Jacobian)%*%K%*%Jacobian
}

    return(max(c(0, sum(diag(solve(nearPD(hess)$mat, K))))))
#}else{  
#    if(length(nonzero.idx)==1) return(0) ## only one beta is nonzero, but pi0=0 too
#    
#    thetas.star=thetas[nonzero.idx[-1],drop=FALSE]
#    h0.i.star=b.i[,nonzero.idx[1]]+h0.i
#    b.i.star=b.i[,nonzero.idx[-1],drop=FALSE]+h0.i-h0.i.star
#    Lik.star=drop(h0.i.star+b.i.star%*%thetas.star)
#    Omega.theta.star=drop(Omega[nonzero.idx,nonzero.idx]%*%thetas[nonzero.idx])
#    C=cbind(-1,diag(1,length(thetas.star)))
#    grads.star=-b.i.star/Lik.star+ matrix(lambda*drop(C%*%Omega.theta.star)/G, G ,length(thetas.star), byrow=TRUE)
#    hess.star=crossprod(b.i.star/Lik.star)+lambda*C%*%Omega[nonzero.idx,nonzero.idx]%*%t(C)
#    K.star=crossprod(grads.star)
#
#    return(max(c(0, sum(diag(solve(nearPD(hess.star)$mat, K.star))))))
#}
}


####################  optimization using constrOptim
#UI=rbind(diag(1,K),-1   )      
#CI=c(rep(0,each=K),-1     )    
#strt=runif(K); strt=strt/sum(strt)*.4
#
#co.fit=constrOptim(strt,NPLL,grad.NPLL,UI,CI)
#grad.NPLL(co.fit$par)
#co.fit$par
#plot(c(1-sum(co.fit$par), co.fit$par/sum(co.fit$par)))
#
#
#co.fit=constrOptim(strt,NPLL,grad.NPLL,UI,CI)
#co.fit=constrOptim(strt,NPLL,NULL,UI,CI)
#curve(d.ncp(x, co.fit$par/sum(co.fit$par)), -5,5,100,add=TRUE)


dev.new(); 
strt=rep(1/K,K); strt=strt/sum(strt)*.5
#strt=runif(K); strt=strt/sum(strt)*.5
sqp.fit=sqp(strt,1e-10)
#curve(d.ncp(x, sqp.fit/sum(sqp.fit)), min(Fstat),max(Fstat),100,add=FALSE,col=1,lwd=1)
hist(ncps,br=50, border=1, xlim=c(min(Fstat)-6*sigs[1],max(Fstat)+6*sigs[1]),pr=TRUE); rug(mus)
#hist(Fstat,pr=TRUE,br=50,sub=lambda); rug(mus)


nics=enps=nic.sd=pi0s=numeric(20)
lambdas=10^seq(-0,2,length=length(nics))
for(i in rev(seq(along=nics))){
    lambda=lambdas[i]
    sqp.fit=sqp(sqp.fit,1e-6)
    enps[i]=enp(sqp.fit)
    tmp=NPLL(sqp.fit,FALSE)
    nics[i]=sum(tmp)+enps[i]
    nic.sd[i]=sd(drop(tmp))*sqrt(G)
    pi0s[i]=1-sum(sqp.fit)
    hist(ncps,pr=TRUE,xlim=c(min(Fstat),max(Fstat)),br=40,sub=lambda,main=enps[i]); rug(mus); 
    curve(d.ncp(x, sqp.fit/sum(sqp.fit)), min(Fstat),max(Fstat),100,add=TRUE,col=4,lwd=2)

#    hist(Fstat,pr=TRUE,br=150,sub=lambda,main=paste('pi0=',1-sum(sqp.fit))); rug(mus); 
#    lines(Fstat[Fstat.ord],d.t(sqp.fit)[Fstat.ord],col=2,lwd=2)
}
plot(log10(lambdas), enps); 
plot(log10(lambdas), pi0s); 
plot(log10(lambdas), ((nics)))
for(i in seq(along=nics))lines(rep(log10(lambdas)[i],2), nics[i]+c(-1,1)*nic.sd[i], lwd=3)
abline(v=log10(lambdas)[which.min(nics)], col=2);
abline(h=nics[which.min(nics)]+c(-1,1)*nic.sd[which.min(nics)], col=4)
abline(v=log10(lambdas)[which.max(nics[nics<=nics[which.min(nics)]+nic.sd[which.min(nics)]])], col=2,lty=2)

lambda=lambdas[which.min(nics)]
#lambda=(lambdas)[which.max(nics[nics<=nics[which.min(nics)]+nic.sd[which.min(nics)]])]
sqp.fit=sqp(sqp.fit,1e-10, 1e10)
hist(ncps,pr=TRUE,xlim=c(min(Fstat),max(Fstat)),br=40,sub=lambda,main=enp(sqp.fit)); rug(mus)
#curve(d.ncp(x, sqp.fit/sum(sqp.fit)), min(Fstat),max(Fstat),100,add=TRUE,col=4,lwd=2)
curve(d.ncp(x, sqp.fit/(sum(sqp.fit))), min(Fstat),max(Fstat),100,add=TRUE,col=4,lwd=2)
#curve(dnorm(x,2), min(Fstat),max(Fstat),100,add=TRUE,col=1,lwd=1)

hist(Fstat,pr=TRUE,br=50,sub=lambda)
lines(Fstat[Fstat.ord],d.t(sqp.fit)[Fstat.ord],col=2,lwd=2)


}