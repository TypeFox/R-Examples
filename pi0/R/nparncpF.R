nparncpF=function(Fstat, df1,df2, ...)
{
    method='SQP'
    if (method=='SQP') {
        nparncpF.sqp(Fstat, df1,df2, ...)
    }
}
nparncpF.sqp = function (Fstat, df1,df2, penalty=3, lambdas=10^seq(-1,8,by=1), starts, IC=c('BIC','CAIC','HQIC','AIC'),
                        K=20, bounds=quantile(Fstat,c(.995)), solver=c('solve.QP','lsei','ipop','LowRankQP'),plotit=FALSE, verbose=FALSE, approx.hess=TRUE, ... )
{
#   source("int.nct.R"); source("laplace.nct.R"); source("saddlepoint.nct.R"); source("dFsnc.mix.R")
    solver=match.arg(solver)
    #penalty=match.arg(penalty)
    solver.package=switch(solver, solve.QP='limSolve', ipop='kernlab', lsei='limSolve',LowRankQP='LowRankQP'
    )
    #library(solver.package)
    #library("Matrix")
    if (K<=0 || length(K)!=1) stop("K should be a positive integer")
    if (any(lambdas<0)) stop("lambdas should be a vector of positive numbers")
    lambdas=sort(lambdas)
#    if (length(bounds)==1) bounds=c(-abs(bounds), abs(bounds))
#    if (diff(bounds)[1]<=0) bounds=bounds[2:1]
    bounds=abs(bounds[1])

	if(length(df1)>1) warning('df1 needs to be a scalar. The first element is used.')
    G=max(c(length(Fstat),length(df2)))
    tmp=seq(0, sqrt(bounds), length=K+1L); 
	mus=c(-rev(tmp[-1]),0,tmp[-1])
	mu2s=tmp^2
	gam2=mu2s[2]

    IC=toupper(IC); IC=match.arg(IC)
    IC.fact=switch(IC, AIC=2, BIC=log(G), CAIC=log(G)+2, HQIC=2*log(log(G)))

    h0.i=df(Fstat,df1,df2) 
    b.i=matrix(0,G,K+1)
    for(k in 1+0:K){
        b.i[,k]=dFsnc.mix(Fstat,df1,df2,mu2s[k], gam2,...)
        if(any(b.i[,k]<0) || any(is.na(b.i[,k]))) {
            warning("Noncentral density unreliable. I switched to exact density function")
            b.i[,k]=dFsnc.mix(Fstat,df1,df2,mu2s[k], gam2,approximation='none')
        }
    }
    if(approx.hess==1) {
        b.approx=as(b.i*(b.i>mean(b.i)), 'dgCMatrix')
    }else if(approx.hess>0 && approx.hess<1) {
        b.approx=as(b.i*(b.i>quantile(b.i,approx.hess)),'dgCMatrix')
    }
    b.i=b.i-h0.i

	symChi=function(k, x){ ## scalar x; vector k
		ans=dchisq(x*x/gam2, df1, mu2s[abs(k)+1L]/gam2)*abs(x)/gam2
		if(df1==1){ 
			idx=(x==0)
			ans[idx]=dnorm(sqrt(mu2s[abs(k[idx])+1L]), , sqrt(gam2))
		}
		return(ans)
	}
    q.mu=outer(-K:K, mus, symChi)
    if(FALSE){
		if(substr(penalty,1,4)=='sqrt') q.mu=sqrt(q.mu)
		if(penalty=='1st.deriv' || penalty=='sqrt.1st.deriv') { ### actually differencing
				diffmat=diag(1,K+2,K+1); diag(diffmat[-1,])=-1
		}else if(penalty=='2nd.deriv' || penalty=='sqrt.2nd.deriv'){
				diffmat=diag(1,K+2,K); diag(diffmat[-1,])=-2; diag(diffmat[-(1:2),])=1
		}else if(penalty=='3rd.deriv' || penalty=='sqrt.3rd.deriv'){
				diffmat=diag(1,K+2,K+1); diag(diffmat[-1,])=-3; diag(diffmat[-(1:2),])=3; diag(diffmat[-(1:3),])=-1
		}
		Omega=tcrossprod(q.mu%*%diffmat)
	}else{
		## FIXME: check penalty
		penalty=as.integer(penalty)
		Omega1=crossprod(apply(q.mu, 1L, diff, difference=penalty))
		L=diag(1, K+1)[1+abs(-K:K),]
		Omega=crossprod(L, Omega1%*%L)
	}

    NPLL=function(thetas, take.sum=TRUE ){ ### depends on h0.i, b.i, lambda, Omega, G
        one_pi0=max(c(sum(thetas), 1e-6))
        Lik=pmax(h0.i+b.i%*%thetas, 1e-8)  #### CHECK negativity PROBLEM
        NLL=-log(Lik)
        if(take.sum) ans=sum(NLL)+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)
        else ans=NLL+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)/G
        ans
    }

    grad.NPLL=function(thetas, take.sum=TRUE) {    ### depends on h0.i, b.i, lambda, Omega, G, K
        one_pi0=max(c(sum(thetas), 1e-6))
        Lik=drop(h0.i+b.i%*%thetas)     #### no log is taken later
        Omega.theta=drop(Omega%*%thetas)
        if(take.sum){
           -drop(crossprod(1/Lik, b.i))+
           lambda/one_pi0/one_pi0*Omega.theta - 
           lambda*drop(crossprod(thetas,Omega.theta))/one_pi0/one_pi0/one_pi0
        }else{  # return a G x k matrix, each row is the gradient based on one data point
            -b.i/Lik+ 
            matrix(lambda/one_pi0/one_pi0*Omega.theta/G, G ,(K+1L), byrow=TRUE)-
            lambda*drop(crossprod(thetas,Omega.theta))/one_pi0/one_pi0/one_pi0/G
        }    
    }
        
    hess.NPLL=function(thetas,approx=FALSE){         ### depends on h0.i, b.i, lambda, Omega, K
        one_pi0=max(c(sum(thetas),1e-6))        ## avoiding division by zero
        Lik= drop(h0.i+b.i%*%thetas)
        Omega.theta=drop(Omega%*%thetas)
        Omega.theta.onep=matrix(Omega.theta,K+1L,K+1L)
        if(approx){
            bl.approx=b.approx/Lik
            h0l=h0.i/Lik
            blh0l=matrix(as(crossprod(bl.approx, h0l),'vector'),K+1,K+1)
            crossprod.term=as.matrix(crossprod(bl.approx))+drop(crossprod(h0l))-blh0l-t(blh0l)
        }else{
            bi.lik=b.i/Lik
            crossprod.term= crossprod(bi.lik)
        }
        crossprod.term + lambda*Omega/one_pi0/one_pi0 - 
            2*lambda*(Omega.theta.onep+t(Omega.theta.onep))/one_pi0/one_pi0/one_pi0+
            3*lambda*matrix(drop(crossprod(thetas,Omega.theta)),K+1L,K+1L)/one_pi0/one_pi0/one_pi0/one_pi0
    }

    C.fctn=function(thetas, eps=1e-8){ ## constraint function s.t. C.fctn(thetas)>=0
        rbind(diag(1,length(thetas)),-1)%*%thetas+rep(c(eps,1-eps),c(length(thetas),1))
    }
    grad.C=function(thetas){  ## grad.C^TRUE thetas + C >=0
        cbind(diag(1,length(thetas)),-1)
    }
    Amat=grad.C(numeric(K+1L)) ## this is the A matrix for quadprog::solve.QP, i.e., t(A)%*%theta>=theta0 linear constraints
                            ## for limSolve::lsei(, this is t(G)


    sqp=function(thetas, conv.f=1e-10, verbose=FALSE, maxiter=1e3) { ## thetas is a starting value
        ## not a general solver; instead, designed for this problem per se
        ## depends on solver, grad.NPLL, hess.NPLL, nearPD, Amat, C.fctn

        npll.last=Inf
        niter=1
        repeat{
            thetas=pmax(thetas,0)
            dvec=-grad.NPLL(thetas)    ## negative gradiant
            fnscale=10^floor(max(log10(abs(dvec))))
            dvec=dvec/fnscale
        #    tmp=drop(1/((1-sum(betas))*F0+F1%*%betas))*(F1-F0)
            Dmat=hess.NPLL(thetas, approx=approx.hess>0)/fnscale
        #    repeat{
        #        eigv=eigen(Dmat,TRUE)
        #        if(tail(eigv$val,1)<1e-5) {
        #            if(verbose)cat('\tmin.eigval=',tail(eigv$val,1),'\trk.D=',qr(Dmat)$rank,fill=TRUE)
        #            Dmat=eigv$vec%*%diag(eigv$val-min(eigv$val)+1.1e-6)%*%t(eigv$vec) 
        #        }else 
        #            break
        #    }
            Dmat=as.matrix(nearPD((Dmat+t(Dmat))/2 )$mat)
        #    Amat=t(UI)
        #    bvec=drop(-UI%*%thetas+CI) -1e-10
            bvec=-C.fctn(thetas)
            qp.sol=if(solver=='solve.QP') {
                        tmpA=try(solve.QP(Dmat,dvec,Amat=Amat,bvec)$solution, silent=TRUE)
                        if(class(tmpA)=='try-error'){ # solver='lsei'
                            tmpA=chol(Dmat); tmpB=.5*drop(forwardsolve(t(tmpA), dvec)); 
                            limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,0,K+1L), F=numeric(0), G=t(Amat), H=bvec, 
                                 Wx=NULL, Wa=NULL, type=1)$X
                        }else tmpA
                   }else if (solver=='lsei') {
                        tmpA=chol(Dmat); tmpB=.5*drop(forwardsolve(t(tmpA), dvec)); 
#                        limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,1,K), F=0, G=t(Amat), H=bvec, 
#                             Wx=NULL, Wa=NULL, type=1)$X
                        limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,0,K+1L), F=numeric(0), G=t(Amat), H=bvec, 
                             Wx=NULL, Wa=NULL, type=1)$X
                   }else if (solver=='ipop') {  ## not working well ## the R translation of the LOQO code is not very honest
                        tmpA=ipop(c=-dvec, H=Dmat, A=t(Amat), b=bvec, l=rep(-1,K+1L), u=rep(1,K+1L), r=rep(1e6,K+1L+1),verb=verbose)
                        if(tmpA@how!='converged'){warning('ipop not converged')}
                        tmpA@primal
                   }else if (solver=='LowRankQP') { ## not working well
                        tmpTransform=solve(tcrossprod(Amat),Amat)
                        Hmat=crossprod(tmpTransform, Dmat%*%tmpTransform)
                        tmpAns=LowRankQP::LowRankQP(Vmat=Hmat,dvec=Hmat%*%bvec-crossprod(tmpTransform,dvec),
                                        Amat=matrix(0,0,K+1),bvec=numeric(0),uvec=rep(1e6,K+1),method='LU')
                        drop(tmpTransform%*%(bvec+tmpAns$alpha))
                   }
            thetas.new=thetas+qp.sol
            repeat{
                if(all(C.fctn(thetas.new)>=0-1e-10) )break  #&& crossprod(Amat,thetas.new)>=bvec
                thetas.new=(thetas+thetas.new)/2
                if(verbose) cat("halving due to constraints",fill=TRUE)
            }
            nhalving=1
            repeat{
                npll.new=NPLL(thetas.new)
                if(npll.new<=npll.last) break
                if(all(thetas==thetas.new)) break
                thetas.new=(thetas+thetas.new)/2
                if(verbose) cat("halving due to increase in objective",fill=TRUE)
                nhalving=nhalving+1
                if(nhalving>50) {warning("step-halving limite reached"); break}
            }
            thetas=thetas.new
            if(verbose)cat('npll.new=', npll.new, '\tpi0=',1-sum(thetas), '\tmin.theta=', min(thetas), '\tfnscale=', fnscale, fill=TRUE)
            if(npll.last-npll.new<conv.f*fnscale) break
            niter=niter+1
            if(niter>maxiter){warning("maximum iteration reached; results from last iteration returned"); break}
            npll.last=npll.new
         }
         if(verbose) cat("DONE!",fill=TRUE)
         return(pmax(thetas,0))
    }

    enp=function(thetas, eps=1e-6)  ## effective number of parameters, ## depends on hess.NPLL, grad.NPLL
    {
        nonzero.idx=which(thetas>eps)

        if(length(nonzero.idx)==0){        # pi0=1
            return(0)
        }
        hess=hess.NPLL(thetas)[nonzero.idx,nonzero.idx,drop=FALSE]
        grads=grad.NPLL(thetas,take.sum=FALSE)[,nonzero.idx,drop=FALSE]
        Kmat=crossprod(grads)

        if( 1-sum(thetas)<eps) {## pi0=0, i.e. sum(thetas)==1; ## this needs reparameterization
            if(length(nonzero.idx)==1) return(0) ## only one beta is nonzero, but pi0=0 too
            Jacobian=rbind(-1,diag(1,length(thetas[nonzero.idx])-1))
            hess=t(Jacobian)%*%hess%*%Jacobian
            Kmat=t(Jacobian)%*%Kmat%*%Jacobian
        }

        return(max(c(0, sum(diag(solve(nearPD(hess)$mat, Kmat))))))
    }

    if(missing(starts)) {
        tmp=parncpF(Fstat,df1,df2,central=FALSE,method='L-BFGS-B')
        starts=rep(coef(tmp)['pi0']/(K+1L),(K+1L))
    }
	if(sum(starts)>1-1e-6) starts=starts/sum(starts)*(1-1e-6)
	starts=pmax(1e-6, starts)

    n.lambda=length(lambdas)
    nics=enps=nic.sd=pi0s=numeric(n.lambda)
    sqp.fit=matrix(NA_real_, n.lambda, (K+1L))
    for(i in n.lambda:1) {
        lambda=lambdas[i]
        sqp.fit[i,]=sqp(if(i==n.lambda) starts else sqp.fit[i+1,], conv.f=1e-6, verbose=verbose)
        enps[i]=enp(sqp.fit[i,])
        tmp=NPLL(sqp.fit[i,],FALSE)
        nics[i]=2*(sum(tmp))+IC.fact*enps[i]
        nic.sd[i]=sd(drop(tmp))*sqrt(G)
        pi0s[i]=1-sum(sqp.fit[i,])
    }

#    if(smooth.enp) {
#          warning("smoothing snp is not well tested")
#library("monoProc")
#        loe=loess(enps~log10(lambdas))
#        mon=mono.1d(list(log10(lambdas), fitted(loe)), bw.nrd0(fitted(loe))/3,mono1='decreasing')
#        enps.smooth=mon@y
#        nics=nics-enps+enps.smooth
#    }

    i.final=which.min(nics); # i.1se=tail(which(nics<=nics[i.final]+nic.sd[i.final] & 1:n.lambda>=i.final),1) # one se rule
    lambda.final=lambdas[i.final]
    beta.final=sqp.fit[i.final,]/sum(sqp.fit[i.final,])


    ll=-NPLL(sqp.fit[i.final,])
    attr(ll, 'df')=enps[i.final]
	attr(ll, 'nobs')=G
    class(ll)='logLik'

	mu.ncp.k=(df1*gam2+mu2s)
	mu.ncp=sum(beta.final * mu.ncp.k)
	mu.ncp2=sum(beta.final*(mu.ncp.k^2+2*gam2*(mu.ncp.k+mu2s)))
    ans=list(pi0=pi0s[i.final], 
		mu.ncp=mu.ncp, sd.ncp= sqrt(mu.ncp2-mu.ncp^2), 
             logLik=ll, enp=enps[i.final], par=sqp.fit[i.final,],lambda=lambdas[i.final],
             gradiant=grad.NPLL(sqp.fit[i.final,]), hessian=hess.NPLL(sqp.fit[i.final,]),

             beta=beta.final, IC=IC, 
             all.mus=mus, gam2=gam2, data=list(Fstat=Fstat, df1=df1, df2=df2), i.final=i.final, all.pi0s=pi0s,
             all.enps=enps, all.thetas=sqp.fit, all.nics=nics, all.nic.sd=nic.sd, all.lambdas=lambdas,nobs=G)
    class(ans)=c('nparncpF','ncpest')
    if(plotit) plot.nparncpF(ans)
    ans
}

fitted.nparncpF=#fitted.values.nparncpF=
function(object, ...)
{
    mu2s=sort(unique(object$all.mus^2))
    nonnull.mat=matrix(NA_real_, length(object$data$Fstat), length(mu2s))
    for(k in 1:length(object$all.gam2s)) nonnull.mat[,k]=dFsnc.mix(object$data$Fstat, object$data$df1,object$data$df2, mu2s,object$gam2,FALSE,...)
    object$pi0*df(object$data$Fstat, object$data$df1,object$data$df2)+(1-object$pi0)*drop(nonnull.mat%*%object$beta)
}

if(FALSE){
summary.nparncpF=function(object,...)
{
    print.nparncpF(object,...)
}
print.nparncpF=function(x,...)
{
    cat('pi0=',x$pi0, fill=TRUE)
    cat('mu.ncp=', x$mu.ncp, fill=TRUE)
    cat('sd.ncp=', x$sd.ncp, fill=TRUE)
    cat("enp=", x$enp, fill=TRUE)
    cat("lambda=", x$lambda, fill=TRUE)
    invisible(x)
}
}

plot.nparncpF=function(x,...)
{
#    dev.new(width=7,heigh=7)
    op=par(mfrow=c(2,2))
#    attach(x)
    n.lambda=length(x$all.lambdas)
    i.2se=which(x$all.nics<=x$all.nics[x$i.final]+x$all.nic.sd[x$i.final]*2)
    if(length(i.2se)==0) i.2se=x$i.final+c(-1,1)
	plot(log10(x$all.lambdas), (x$all.nics), 
        ylim=range(c(x$all.nics+x$all.nic.sd, x$all.nics-x$all.nic.sd)[c(i.2se,n.lambda+i.2se)]), 
        xlab='log10(lambda)', ylab='NIC')
    for(i in 1:n.lambda) lines(rep(log10(x$all.lambdas)[i],2), x$all.nics[i]+c(-1,1)*x$all.nic.sd[i], lwd=3)
    abline(v=log10(x$all.lambdas[x$i.final]), col=2);
    abline(h=x$all.nics[x$i.final]+c(-1,1)*x$all.nic.sd[x$i.final], col=4)
#    abline(v=log10(x$all.lambdas[i.1se]), col=2,lty=2)
	title(sub=paste('lambda =',round(x$lambda,3),sep=''))
	
#    plot(log10(x$all.lambdas), x$all.enps); abline(v=log10(x$all.lambdas[c(x$i.final,i.1se)]), xlab='log10(lambda)', ylab='ENP', lty=1:2)
#    plot(log10(x$all.lambdas), x$all.pi0s); abline(v=log10(x$all.lambdas[c(x$i.final,i.1se)]), xlab='log10(lambda)', ylab='pi0', lty=1:2)
    plot(log10(x$all.lambdas), x$all.enps,xlab='log10(lambda)',ylab='effective # parameters', ylim=c(0,max(x$all.enps[i.2se]))); 
        abline(v=log10(x$all.lambdas[c(x$i.final)]), xlab='log10(lambda)', ylab='ENP', lty=1)
    title(sub=paste('df =', round(x$enp,3),sep=''))
	
	plot(log10(x$all.lambdas), x$all.pi0s, xlab='log10(lambda)', ylab='pi0',ylim=0:1); 
        abline(v=log10(x$all.lambdas[c(x$i.final)]), xlab='log10(lambda)', ylab='pi0', lty=1)
	title(sub=paste('pi0 =', round(x$pi0,3),sep=''))
	
	mu2s=sort(unique(x$all.mus^2))
    d.ncp=function(xx)        # p in the paper
    {   ## depends on mu2s, x
		schisq.x=outer(xx/x$gam2, mu2s/x$gam2, dchisq, df=x$data$df1)/x$gam2
		ans=drop(schisq.x %*% x$beta)
		if(any(xx==0) && x$data$df1<2){
			ans[xx==0]=Inf
		}
		ans
    }
	curve(d.ncp, 0, max(mu2s)*1.1, col=4,lwd=2, xlab=expression(delta),ylab='density')
	lines(c(-max(mu2s)*.05,0),c(0,0), col=4, lwd=2)
#    detach(x)
    rug(mu2s)
    title(sub=paste("mean =",round(x$mu.ncp,3),"; sd =",round(x$sd.ncp,3),sep=''))
    par(op)
    invisible(x)
}

if(FALSE) {
	#library(pi0)
	environment(nparncpF)=environment(nparncpF.sqp)=environment(nparncpt)
	set.seed(243)
	ncpt=c(rep(0,.31e3),rnorm(.69e3, .5))
	df1=1; df2=4
	tstat=rt(1e4, df2, ncpt)
	fstat=tstat^2
	(npf=nparncpF.sqp(fstat,df1,df2)); plot(npf)
	(npt=nparncpt(tstat, df2, plotit=TRUE))

	if(FALSE){## old way of defining bases
		nparncpF.sqp = function (Fstat, df1,df2, penalty=c('sqrt.3rd.deriv','sqrt.2nd.deriv','sqrt.1st.deriv','3rd.deriv','2nd.deriv','1st.deriv'), lambdas=10^seq(-1,8,by=1), starts, IC=c('BIC','CAIC','HQIC','AIC'),
                        K=20, bounds=quantile(Fstat,c(.99)), solver=c('solve.QP','lsei','ipop','LowRankQP'),plotit=FALSE, verbose=FALSE, approx.hess=TRUE, ... )
{
#   source("int.nct.R"); source("laplace.nct.R"); source("saddlepoint.nct.R"); source("dFsnc.mix.R")
    solver=match.arg(solver)
    penalty=match.arg(penalty)
    solver.package=switch(solver, solve.QP='limSolve', ipop='kernlab', lsei='limSolve',LowRankQP='LowRankQP'
    )
    #library(solver.package)
    #library("Matrix")
    if (K<=0 || length(K)!=1) stop("K should be a positive integer")
    if (any(lambdas<0)) stop("lambdas should be a vector of positive numbers")
    lambdas=sort(lambdas)
#    if (length(bounds)==1) bounds=c(-abs(bounds), abs(bounds))
#    if (diff(bounds)[1]<=0) bounds=bounds[2:1]
    bounds=abs(bounds[1])

    G=max(c(length(Fstat),length(df1),length(df2)))
    mus=seq(0, bounds, length=K+2)
    mus[1]=mus[2]/10
    gam2s=numeric(K)
    for(k in 1:K) gam2s[k]=optimize(function(g2)abs(g2*qchisq(.95,df1,mus[k]/g2-df1)-mus[k+2]),c(1e-3,mus[k]/df1))$min
#                get.gam2=function(){
#                  gam2s=numeric(K)
#                  for(k in 1:K){
#                    roots=polyroot(c(mus[2]^2,-4*mus[k],2*df1))
#                    Re.rt=Re(roots); Im.rt=Im(roots)=c(0,0)
#                    nroots=sum(Im.rt==0)
#                    if(nroots==2){
#                        good=Re.rt>0 & Re.rt<=mus[k]/df1
#                        if(sum(good)==1){ gam2s[k]=Re.rt[good]; next}
#                        if(sum(good)==2) browser()  ## this is unreasonable result
#                    }else if(nroots==1){
#                        Re.rt=Re.rt[Im.rt==0]
#                        if( Re.rt>0 & Re.rt<=mus[k]/df1 ) {gam2s[k]=Re.rt; next}
#                    }
#                    cat(k, nroots, fill=T)
#                    gam2s[k]=optimize(function(g2)abs(g2*qchisq(.95,df1,mus[k]/g2-df1)-mus[k+2]),c(1e-3,mus[k]/df1))$min
#                  }
#                  gam2s
#                }
#                gam2s=get.gam2()

    IC=toupper(IC); IC=match.arg(IC)
    IC.fact=switch(IC, AIC=2, BIC=log(G), CAIC=log(G)+2, HQIC=2*log(log(G)))

    h0.i=df(Fstat,df1,df2) 
    b.i=matrix(0,G,K)
    for(k in 1:K){
        b.i[,k]=dFsnc.mix(Fstat,df1,df2,mus[k], gam2s[k],...)
        if(any(b.i[,k]<0) || any(is.na(b.i[,k]))) {
            warning("Noncentral density unreliable. I switched to exact density function")
            b.i[,k]=dFsnc.mix(Fstat,df1,df2,mus[k], gam2s[k],approximation='none')
        }
    }
    if(approx.hess==1) {
        b.approx=as(b.i*(b.i>mean(b.i)), 'dgCMatrix')
    }else if(approx.hess>0 && approx.hess<1) {
        b.approx=as(b.i*(b.i>quantile(b.i,approx.hess)),'dgCMatrix')
    }
    b.i=b.i-h0.i


    q.mu=outer(1:K, mus, function(k, u) dchisq(u/gam2s[k], df1, mus[k]/gam2s[k]-df1)/gam2s[k])
    if(substr(penalty,1,4)=='sqrt') q.mu=sqrt(q.mu)
    if(penalty=='1st.deriv' || penalty=='sqrt.1st.deriv') { ### actually differencing
            diffmat=diag(1,K+2,K+1); diag(diffmat[-1,])=-1
    }else if(penalty=='2nd.deriv' || penalty=='sqrt.2nd.deriv'){
            diffmat=diag(1,K+2,K); diag(diffmat[-1,])=-2; diag(diffmat[-(1:2),])=1
    }else if(penalty=='3rd.deriv' || penalty=='sqrt.3rd.deriv'){
            diffmat=diag(1,K+2,K+1); diag(diffmat[-1,])=-3; diag(diffmat[-(1:2),])=3; diag(diffmat[-(1:3),])=-1
    }
    Omega=tcrossprod(q.mu%*%diffmat)

    NPLL=function(thetas, take.sum=TRUE ){ ### depends on h0.i, b.i, lambda, Omega, G
        one_pi0=max(c(sum(thetas), 1e-6))
        Lik=pmax(h0.i+b.i%*%thetas, 1e-8)  #### CHECK negativity PROBLEM
        NLL=-log(Lik)
        if(take.sum) ans=sum(NLL)+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)
        else ans=NLL+.5*lambda*drop(thetas%*%Omega%*%thetas)/(one_pi0*one_pi0)/G
        ans
    }

    grad.NPLL=function(thetas, take.sum=TRUE) {    ### depends on h0.i, b.i, lambda, Omega, G, K
        one_pi0=max(c(sum(thetas), 1e-6))
        Lik=drop(h0.i+b.i%*%thetas)     #### no log is taken later
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
        
    hess.NPLL=function(thetas,approx=FALSE){         ### depends on h0.i, b.i, lambda, Omega, K
        one_pi0=max(c(sum(thetas),1e-6))        ## avoiding division by zero
        Lik= drop(h0.i+b.i%*%thetas)
        Omega.theta=drop(Omega%*%thetas)
        Omega.theta.onep=matrix(Omega.theta,K,K)
        if(approx){
            bl.approx=b.approx/Lik
            h0l=h0.i/Lik
            blh0l=matrix(as(crossprod(bl.approx, h0l),'vector'),K,K)
            crossprod.term=as.matrix(crossprod(bl.approx))+drop(crossprod(h0l))-blh0l-t(blh0l)
        }else{
            bi.lik=b.i/Lik
            crossprod.term= crossprod(bi.lik)
        }
        crossprod.term + lambda*Omega/one_pi0/one_pi0 - 
            2*lambda*(Omega.theta.onep+t(Omega.theta.onep))/one_pi0/one_pi0/one_pi0+
            3*lambda*matrix(drop(crossprod(thetas,Omega.theta)),K,K)/one_pi0/one_pi0/one_pi0/one_pi0
    }

    C.fctn=function(thetas, eps=1e-8){ ## constraint function s.t. C.fctn(thetas)>=0
        rbind(diag(1,length(thetas)),-1)%*%thetas+rep(c(eps,1-eps),c(length(thetas),1))
    }
    grad.C=function(thetas){  ## grad.C^TRUE thetas + C >=0
        cbind(diag(1,length(thetas)),-1)
    }
    Amat=grad.C(numeric(K)) ## this is the A matrix for quadprog::solve.QP, i.e., t(A)%*%theta>=theta0 linear constraints
                            ## for limSolve::lsei(, this is t(G)


    sqp=function(thetas, conv.f=1e-10, verbose=FALSE, maxiter=1e3) { ## thetas is a starting value
        ## not a general solver; instead, designed for this problem per se
        ## depends on solver, grad.NPLL, hess.NPLL, nearPD, Amat, C.fctn

        npll.last=Inf
        niter=1
        repeat{
            thetas=pmax(thetas,0)
            dvec=-grad.NPLL(thetas)    ## negative gradiant
            fnscale=10^floor(max(log10(abs(dvec))))
            dvec=dvec/fnscale
        #    tmp=drop(1/((1-sum(betas))*F0+F1%*%betas))*(F1-F0)
            Dmat=hess.NPLL(thetas, approx=approx.hess>0)/fnscale
        #    repeat{
        #        eigv=eigen(Dmat,TRUE)
        #        if(tail(eigv$val,1)<1e-5) {
        #            if(verbose)cat('\tmin.eigval=',tail(eigv$val,1),'\trk.D=',qr(Dmat)$rank,fill=TRUE)
        #            Dmat=eigv$vec%*%diag(eigv$val-min(eigv$val)+1.1e-6)%*%t(eigv$vec) 
        #        }else 
        #            break
        #    }
            Dmat=as.matrix(nearPD((Dmat+t(Dmat))/2 )$mat)
        #    Amat=t(UI)
        #    bvec=drop(-UI%*%thetas+CI) -1e-10
            bvec=-C.fctn(thetas)
            qp.sol=if(solver=='solve.QP') {
                        tmpA=try(solve.QP(Dmat,dvec,Amat=Amat,bvec)$solution, silent=TRUE)
                        if(class(tmpA)=='try-error'){ # solver='lsei'
                            tmpA=chol(Dmat); tmpB=.5*drop(forwardsolve(t(tmpA), dvec)); 
                            limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,0,K), F=numeric(0), G=t(Amat), H=bvec, 
                                 Wx=NULL, Wa=NULL, type=1)$X
                        }else tmpA
                   }else if (solver=='lsei') {
                        tmpA=chol(Dmat); tmpB=.5*drop(forwardsolve(t(tmpA), dvec)); 
#                        limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,1,K), F=0, G=t(Amat), H=bvec, 
#                             Wx=NULL, Wa=NULL, type=1)$X
                        limSolve::lsei(A=tmpA, B=tmpB, E=matrix(0,0,K), F=numeric(0), G=t(Amat), H=bvec, 
                             Wx=NULL, Wa=NULL, type=1)$X
                   }else if (solver=='ipop') {  ## not working well ## the R translation of the LOQO code is not very honest
                        tmpA=ipop(c=-dvec, H=Dmat, A=t(Amat), b=bvec, l=rep(-1,K), u=rep(1,K), r=rep(1e6,K+1),verb=verbose)
                        if(tmpA@how!='converged'){warning('ipop not converged')}
                        tmpA@primal
                   }else if (solver=='LowRankQP') { ## not working well
                        tmpTransform=solve(tcrossprod(Amat),Amat)
                        Hmat=crossprod(tmpTransform, Dmat%*%tmpTransform)
                        tmpAns=LowRankQP::LowRankQP(Vmat=Hmat,dvec=Hmat%*%bvec-crossprod(tmpTransform,dvec),
                                        Amat=matrix(0,0,K+1),bvec=numeric(0),uvec=rep(1e6,K+1),method='LU')
                        drop(tmpTransform%*%(bvec+tmpAns$alpha))
                   }
            thetas.new=thetas+qp.sol
            repeat{
                if(all(C.fctn(thetas.new)>=0-1e-10) )break  #&& crossprod(Amat,thetas.new)>=bvec
                thetas.new=(thetas+thetas.new)/2
                if(verbose) cat("halving due to constraints",fill=TRUE)
            }
            nhalving=1
            repeat{
                npll.new=NPLL(thetas.new)
                if(npll.new<=npll.last) break
                if(all(thetas==thetas.new)) break
                thetas.new=(thetas+thetas.new)/2
                if(verbose) cat("halving due to increase in objective",fill=TRUE)
                nhalving=nhalving+1
                if(nhalving>50) {warning("step-halving limite reached"); break}
            }
            thetas=thetas.new
            if(verbose)cat('npll.new=', npll.new, '\tpi0=',1-sum(thetas), '\tmin.theta=', min(thetas), '\tfnscale=', fnscale, fill=TRUE)
            if(npll.last-npll.new<conv.f*fnscale) break
            niter=niter+1
            if(niter>maxiter){warning("maximum iteration reached; results from last iteration returned"); break}
            npll.last=npll.new
         }
         if(verbose) cat("DONE!",fill=TRUE)
         return(pmax(thetas,0))
    }

    enp=function(thetas, eps=1e-6)  ## effective number of parameters, ## depends on hess.NPLL, grad.NPLL
    {
        nonzero.idx=which(thetas>eps)

        if(length(nonzero.idx)==0){        # pi0=1
            return(0)
        }
        hess=hess.NPLL(thetas)[nonzero.idx,nonzero.idx,drop=FALSE]
        grads=grad.NPLL(thetas,take.sum=FALSE)[,nonzero.idx,drop=FALSE]
        Kmat=crossprod(grads)

        if( 1-sum(thetas)<eps) {## pi0=0, i.e. sum(thetas)==1; ## this needs reparameterization
            if(length(nonzero.idx)==1) return(0) ## only one beta is nonzero, but pi0=0 too
            Jacobian=rbind(-1,diag(1,length(thetas[nonzero.idx])-1))
            hess=t(Jacobian)%*%hess%*%Jacobian
            Kmat=t(Jacobian)%*%Kmat%*%Jacobian
        }

        return(max(c(0, sum(diag(solve(nearPD(hess)$mat, Kmat))))))
    }

    if(missing(starts)) {
        tmp=parncpF(Fstat,df1,df2,central=FALSE)
        starts=rep(coef(tmp)['pi0']/K,K)
        if(sum(starts)>1-1e-6) starts=starts/sum(starts)*(1-1e-6)
    }

    n.lambda=length(lambdas)
    nics=enps=nic.sd=pi0s=numeric(n.lambda)
    sqp.fit=matrix(NA_real_, n.lambda, K)
    for(i in n.lambda:1) {
        lambda=lambdas[i]
        sqp.fit[i,]=sqp(if(i==n.lambda) starts else sqp.fit[i+1,], conv.f=1e-6, verbose=verbose)
        enps[i]=enp(sqp.fit[i,])
        tmp=NPLL(sqp.fit[i,],FALSE)
        nics[i]=2*(sum(tmp))+IC.fact*enps[i]
        nic.sd[i]=sd(drop(tmp))*sqrt(G)
        pi0s[i]=1-sum(sqp.fit[i,])
    }

#    if(smooth.enp) {
#          warning("smoothing snp is not well tested")
#library("monoProc")
#        loe=loess(enps~log10(lambdas))
#        mon=mono.1d(list(log10(lambdas), fitted(loe)), bw.nrd0(fitted(loe))/3,mono1='decreasing')
#        enps.smooth=mon@y
#        nics=nics-enps+enps.smooth
#    }

    i.final=which.min(nics); # i.1se=tail(which(nics<=nics[i.final]+nic.sd[i.final] & 1:n.lambda>=i.final),1) # one se rule
    lambda.final=lambdas[i.final]
    beta.final=sqp.fit[i.final,]/sum(sqp.fit[i.final,])


    ll=-NPLL(sqp.fit[i.final,])
    attr(ll, 'df')=enps[i.final]
	attr(ll, 'nobs')=G
    class(ll)='logLik'

    ans=list(pi0=pi0s[i.final], mu.ncp=beta.final%*%mus[1:K], sd.ncp= sqrt(beta.final%*%(mus[1:K]^2+2*gam2s^2*(2*mus[1:K]/gam2s-df1))-(beta.final%*%mus[1:K])^2), 
             logLik=ll, enp=enps[i.final], par=sqp.fit[i.final,],lambda=lambdas[i.final],
             gradiant=grad.NPLL(sqp.fit[i.final,]), hessian=hess.NPLL(sqp.fit[i.final,]),

             beta=beta.final, IC=IC, 
             all.mus=mus, all.gam2s=gam2s, data=list(Fstat=Fstat, df1=df1, df2=df2), i.final=i.final, all.pi0s=pi0s,
             all.enps=enps, all.thetas=sqp.fit, all.nics=nics, all.nic.sd=nic.sd, all.lambdas=lambdas,nobs=G)
    class(ans)=c('nparncpF','ncpest')
    if(plotit) plot.nparncpF(ans)
    ans
	}

	plot.nparncpF=function(x,...)
{
#    dev.new(width=7,heigh=7)
    op=par(mfrow=c(2,2))
#    attach(x)
    n.lambda=length(x$all.lambdas)
    i.1se=tail(which(x$all.nics<=x$all.nics[x$i.final]+x$all.nic.sd[x$i.final] & 1:n.lambda>=x$i.final),1) ## one se rule
    plot(log10(x$all.lambdas), (x$all.nics), 
        ylim=range(c(x$all.nics+x$all.nic.sd, x$all.nics-x$all.nic.sd)[c(1:min(i.1se,x$i.final+1,n.lambda),n.lambda+1:min(i.1se,x$i.final+1,n.lambda))]), 
        xlab='log10(lambda)', ylab='NIC')
    for(i in 1:n.lambda) lines(rep(log10(x$all.lambdas)[i],2), x$all.nics[i]+c(-1,1)*x$all.nic.sd[i], lwd=3)
    abline(v=log10(x$all.lambdas[x$i.final]), col=2);
    abline(h=x$all.nics[x$i.final]+c(-1,1)*x$all.nic.sd[x$i.final], col=4)
#    abline(v=log10(x$all.lambdas[i.1se]), col=2,lty=2)

#    plot(log10(x$all.lambdas), x$all.enps); abline(v=log10(x$all.lambdas[c(x$i.final,i.1se)]), xlab='log10(lambda)', ylab='ENP', lty=1:2)
#    plot(log10(x$all.lambdas), x$all.pi0s); abline(v=log10(x$all.lambdas[c(x$i.final,i.1se)]), xlab='log10(lambda)', ylab='pi0', lty=1:2)
    plot(log10(x$all.lambdas), x$all.enps,xlab='log10(lambda)',ylab='effective # parameters'); 
        abline(v=log10(x$all.lambdas[c(x$i.final)]), xlab='log10(lambda)', ylab='ENP', lty=1)
    plot(log10(x$all.lambdas), x$all.pi0s, xlab='log10(lambda)', ylab='pi0'); 
        abline(v=log10(x$all.lambdas[c(x$i.final)]), xlab='log10(lambda)', ylab='pi0', lty=1)

    d.ncp=function(xx)        # p in the paper
    {   ## depends on mus, sigs
        z=function(k, u) dchisq(u/x$all.gam2s[k], x$data$df1, x$all.mus[k]/x$all.gam2s[k]-x$data$df1)/x$all.gam2s[k]
        qx=outer(1:(length(x$all.mus)-2), xx,  z)
        drop(x$beta %*% qx)
    }
    curve(d.ncp, min(x$all.mus),max(x$all.mus),100,col=4,lwd=2, xlab=expression(delta),ylab='density')
#    detach(x)
    rug(x$all.mus)
    par(op)
    invisible(x)
}
	fitted.nparncpF=#fitted.values.nparncpF=
function(object, ...)
{
    
    nonnull.mat=matrix(NA_real_, length(object$data$Fstat), length(object$all.mus)-2)
    for(k in 1:length(object$all.gam2s)) nonnull.mat[,k]=dFsnc.mix(object$data$Fstat, object$data$df1,object$data$df2, object$all.mus[k],object$all.gam2s[k],FALSE,...)
    object$pi0*df(object$data$Fstat, object$data$df1,object$data$df2)+(1-object$pi0)*drop(nonnull.mat%*%object$beta)
}
}

}

