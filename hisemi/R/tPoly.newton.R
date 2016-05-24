########## EM algorithm followed by Newton-type optimization
tPoly.newton=function(tstat,x,df,
     #spar=c(10^seq(-1,8,length=30), Inf), nknots=100, 
     starts,
     #tuning.method=c('NIC','CV') ,#'GCV','BIC','CAIC','HQIC'),
     #cv.fold=5, 
     pen.order=1, #poly.degree=pen.order*2-1,
     optim.method=c('nlminb',"BFGS","CG","L-BFGS-B","Nelder-Mead", "SANN", 'NR'),
     #logistic.correction=TRUE,
     #em.iter.max=10, 
     #em.beta.iter.max=1,
     newton.iter.max=1500,
     scale.conv=1e-3, lfdr.conv=1e-3, NPLL.conv=1e-3, 
     debugging=FALSE, plotit=TRUE,...)
{
    if(pen.order==1) return(scaledTMix.null(tstat, df))

    if(debugging)options(error=quote(dump.frames("testdump", TRUE)))  ## this will save the objects when error occurs
    pen.order=round(pen.order)
    poly.degree=round(pen.order)    #round(poly.degree)
    #cv.fold=round(cv.fold)
    
    ##### initializing central t and non-central t
    dt0.all=dt(tstat,df)
    ###### saturated fit           
    logLik.saturated=scaledTMix.sat(tstat,df)

    optim.method=match.arg(optim.method)
    #stopifnot(all(spar>=0))
    #spar=sort(spar)
    #spar=unique(c(spar,Inf))
    spar=Inf

    ##### initializing CV
    #tuning.method=match.arg(tuning.method)
    #if(tuning.method!='CV'){
        tuning.method='NIC'
        cv.fold=1
        noCV=TRUE
    #}else{
    #    cv.fold=ceiling(cv.fold)
    #    if(cv.fold<2){
    #        cv.fold=5
    #        warning(paste("when using tuning.method='CV', cv.fold has to be at least 2. It is set to", cv.fold)
    #        )
    #    }
    #    noCV=FALSE
    #}
    tstat.all=tstat
    G=G.all=length(tstat)
    n.spar=length(spar)
    cv.grp=rep(1:cv.fold,length=G)
    cv.grp=sample(cv.grp)
    criterion.mean=enps=enps0=criterion.var=matrix(, cv.fold, n.spar, 
        dimnames=list(paste("cv.",1:cv.fold,sep=''), paste("log10spar.",log10(spar),sep='')))


    ######## creating quadratic B-spline basis
    #nknots=min(c(nknots, G.all-2))
    if(!inherits(x,'matrix')) x=as.matrix(x)
    n.vars=ncol(x);
    
    #    require(splines)
    #    library(Matrix)
    #    require(fda)

    H.all=Matrix(matrix(1,G.all,1),sparse=TRUE);   ########### intercept term is not penalized 
    j.all=0
    #Pen.mat=matrix(0,1,1)       ########### intercept term is not penalized
    for(i in 1:n.vars){
        #knots.all=quantile(unique(x[,i]), 1:nknots/(nknots+1))
        #H.i=bs(x[,i],knots=knots.all, degree=poly.degree, intercept=FALSE)  ########### intercept term is not penalized
        H.i=poly(x[,i], degree=pen.order-1)
        H.all=cBind(H.all,H.i[,])
        j.all=c(j.all, rep(i, ncol(H.i)))
        ######## derivative penalties 
        #if(pen.order*2-1==poly.degree){
        #  Pen.i=OsplinePen(range(x[,i]), knots.all, pen.order)[-1,-1]
        #}else if (pen.order != poly.degree) {
        #  Hobj.i=create.bspline.basis(breaks=c(min(x[,i]),knots.all,max(x[,i])),norder=poly.degree+1)  ## norder=3 means quadratic
        #  Pen.i=bsplinepen(Hobj.i, Lfdobj=pen.order)[-1,-1] ########### intercept term is not penalized
        #}else {
        #  Hobj.i=create.bspline.basis(breaks=c(min(x[,i]),knots.all,max(x[,i])),norder=poly.degree+1)  ## norder=3 means quadratic
        #  Pen.i=bsplinepen.fda(Hobj.i, Lfdobj=pen.order)[-1,-1] ########### intercept term is not penalized
        #}
        
        #Pen.mat=directSum(Pen.mat,Pen.i)
    }
    #Pen.mat=Matrix(Pen.mat,sparse=TRUE)
    all.parms=matrix(0.0, 1+ncol(H.all), n.spar)

    ############ penalized logLik of data and logLik of test data
    #    logit=make.link("logit")$linkfun  ## this is the logit link function in glm of package stats
    nLogLik.pen=function(parms, with.grad.hess=FALSE, ...) ## also depends on dt0, tstat,  H
    {   a=parms[1]; scale.fact=1+exp(a);  # scale.fact=parms[1]; alternative parameterization to remove boundary
        betas=parms[-1]
        pi0s=1/(1+exp(H%*%betas))
        ll=sum(log(
            pi0s*dt0+(1-pi0s)*dt(tstat/scale.fact,df)/scale.fact
        ))
        ans=-ll#+.5*drop(betas%*%spar.Pen.mat%*%betas)
        if(with.grad.hess){
            dt1=dt(tstat/scale.fact,df)/scale.fact
            f=drop(pi0s*dt0+(1-pi0s)*dt1)
            lfdrs=drop(pi0s*dt0/f)
            w.beta=drop((dt0-dt1)/f*(1-pi0s)*pi0s)
            w.beta.H=w.beta*H
            d.dbetas=colSums(w.beta.H)#+ drop(betas%*%spar.Pen.mat)
            d.dr=sum(-(tstat*tstat-scale.fact*scale.fact)*exp(log(1-1/scale.fact)-log(tstat*tstat/df+scale.fact*scale.fact)+log(1-lfdrs)))
            attr(ans,'gradient')=c(d.dr,       d.dbetas)

            W.beta.beta=w.beta*drop(pi0s*lfdrs-(1-pi0s)*(1-lfdrs))
            Wbb.H=W.beta.beta*H
            HWH=as(crossprod(H,Wbb.H),'symmetricMatrix')
            hess.bb=HWH#+spar.Pen.mat
        
            t2_s2=(tstat*tstat-scale.fact*scale.fact)
            t2vs2=(tstat*tstat/df+scale.fact*scale.fact)

            d.dr.i=drop(-t2_s2*exp(log(1-1/scale.fact)-log(t2vs2)+log(1-lfdrs)))
            w.br=drop(lfdrs*d.dr.i)
            d.dbr=drop(crossprod(H,w.br))

            d2.dr2=sum((scale.fact-1)^2/scale.fact*(1-lfdrs)*t2_s2/t2vs2*(
                    -lfdrs/scale.fact*t2_s2/t2vs2+1/scale.fact/(1-scale.fact)
                    +2*scale.fact/t2_s2 +2*scale.fact/t2vs2
               ))
            attr(ans,'hessian')=as(rBind(c(d2.dr2, drop(d.dbr)),cBind(drop(d.dbr),hess.bb)),'symmetricMatrix')
        }
        ans
    }

    deriv.nLogLik.pen=function(parms, return.K=FALSE, ...)  ## also depends on dt0, tstat,, H
    {   r=parms[1]; scale.fact=1+exp(r);  # scale.fact=parms[1]; alternative parameterization to remove boundary
        betas=parms[-1]
        pi0s=drop(1/(1+exp(H%*%betas)))
        dt1=dt(tstat/scale.fact,df)/scale.fact
        f=drop(pi0s*dt0+(1-pi0s)*dt1)
        lfdrs=drop(pi0s*dt0/f)

        w.beta=drop((dt0-dt1)/f*(1-pi0s)*pi0s)
        w.beta.H=w.beta*H
        d.dbetas=colSums(w.beta.H)#+ drop(betas%*%spar.Pen.mat)
        #        d.dr.i=-exp(-.5*log(df)-lbeta(.5,df/2)+(df-1)*log(scale.fact)+r)*
        #                (exp(
        #                    log(1-pi0s)-log(f)-(df+3)/2*log(tstat*tstat/df+scale.fact*scale.fact)
        #                )*(tstat*tstat-scale.fact*scale.fact))
        d.dr.i=-(tstat*tstat-scale.fact*scale.fact)*exp(log(1-1/scale.fact)-log(tstat*tstat/df+scale.fact*scale.fact)+log(1-lfdrs))
        d.dr=sum(d.dr.i)

        ans=c(d.dr,
              d.dbetas)
#        if(return.K){
#            #               w.beta.H.pen=w.beta.H+matrix(1/length(dt1),length(dt1),1)%*%(betas%*%spar.Pen.mat)   ## this is faster
#            #               dr.w.beta.H.pen=cBind(d.dr.i,    w.beta.H.pen)
#            #               K.full=crossprod(dr.w.beta.H.pen)/length(dt1)    ## -tcrossprod(ans/length(dt1)) ## no need to subtract zero matrix
#            #               attr(ans,'K')=K.full
#                           #### below is equivalent but much faster
#               K.vec=drop(d.dr.i%*%w.beta.H)
#               K=rBind(c(sum(d.dr.i*d.dr.i), (K.vec)), cBind(K.vec, crossprod(w.beta.H)#-crossprod(betas%*%spar.Pen.mat)/G
#               ) ) #/G
#               attr(ans,'K')=as(as(K,'symmetricMatrix'),'sparseMatrix')
#        }
        ans
    }

    hess.nLogLik.pen=function(parms, return.dense=FALSE, ...)  ## also depends on dt0, tstat, , H
    {   r=parms[1]; scale.fact=1+exp(r);  # scale.fact=parms[1]; alternative parameterization to remove boundary
        betas=parms[-1]
        pi0s=drop(1/(1+exp(H%*%betas)))
        dt1=dt(tstat/scale.fact,df)/scale.fact
        f=drop(pi0s*dt0+(1-pi0s)*dt1)
        lfdrs=drop(pi0s*dt0/f)

        w.beta=drop((dt0-dt1)/f*(1-pi0s)*pi0s)
        W.beta.beta=w.beta*drop(pi0s*lfdrs-(1-pi0s)*(1-lfdrs))
        #        W.beta.beta=drop((dt0-dt1)*(f*(2*pi0s-1)+(dt0-dt1)*pi0s*(1-pi0s))*pi0s*(1-pi0s)/f/f)
        Wbb.H=W.beta.beta*H
        HWH=as(crossprod(H,Wbb.H),'symmetricMatrix')
        hess.bb=HWH#+spar.Pen.mat
        hess.bb

        t2_s2=(tstat*tstat-scale.fact*scale.fact)
        t2vs2=(tstat*tstat/df+scale.fact*scale.fact)

        d.dr.i=drop(-t2_s2*exp(log(1-1/scale.fact)-log(t2vs2)+log(1-lfdrs)))
        w.br=drop(lfdrs*d.dr.i)
        d.dbr=drop(crossprod(H,w.br))

        d2.dr2=sum((scale.fact-1)^2/scale.fact*(1-lfdrs)*t2_s2/t2vs2*(
                    -lfdrs/scale.fact*t2_s2/t2vs2+1/scale.fact/(1-scale.fact)
                    +2*scale.fact/t2_s2 +2*scale.fact/t2vs2
               ))
        
        
        ans=rBind(c(d2.dr2, drop(d.dbr)),cBind(drop(d.dbr),hess.bb))
        if(return.dense) as.matrix(ans) else as(ans, 'symmetricMatrix')

    }

    logLiktest=function(betas, scale.fact) ## also depends on H.test, dt0.test, and tstat.test
    {
        fx.test=drop(H.test%*%betas)
        pi0.test=1/(1+exp(fx.test))
        log(
            pi0.test*dt0.test+(1-pi0.test)*dt(tstat.test/scale.fact,df)/scale.fact
        )
    }

    #require(MASS)


    ############ setting staring values
        if(missing(starts)){
            #if(pen.order==1){
            #    null.fit=scaledTMix.null(tstat.all,df)
            #    starts=c(coef(null.fit), rep(0, ncol(H.all)-1))
            #}else {#.NotYetImplemented() ## 
                ########## sequentially increase order of orthogonal polynomials
                tmp.coef=coef( scaledTMix.null(tstat.all,df) )
                #starts=c(coef(null.fit), rep(0, ncol(H.all)-1))
                cur.pen.ord=1
                starts=c(tmp.coef, rep(0,n.vars))
                repeat{
                   cur.pen.ord=cur.pen.ord+1
                   if(cur.pen.ord==pen.order) break

                   tmp = Recall(tstat,x,df,                     
                         starts,
                         pen.order=cur.pen.ord,
                         optim.method,
                         newton.iter.max,
                         scale.conv, lfdr.conv, NPLL.conv, 
                         debugging=FALSE, plotit=FALSE, ...)
                   tmp.coef=coef(tmp)
                   starts=c(tmp.coef[1:2], as.vector(rbind(matrix(tmp.coef[-2:-1], ncol=n.vars), 0)))
                }
            #}
        }


    ############ crossvalidation loop
  for(cv.i in 1:cv.fold){
    cv.idx=which(cv.grp!=cv.i)
    G=length(cv.idx)
    if(length(cv.idx)==0) {cv.idx=1:G.all; G=G.all} ## NO CV at all; using NIC or GCV
    dt0=dt0.all[cv.idx]
    tstat=tstat.all[cv.idx]
    dt0.test=if(noCV) dt0 else dt0.all[-cv.idx]
    tstat.test=if(noCV) tstat else tstat.all[-cv.idx]
    H=H.all[cv.idx,]
    H.test=if(noCV) H else H.all[-cv.idx,]
    #    n2ll.sat.test=attr(logLik.saturated,'n2ll')[if(noCV)cv.idx else -cv.idx]
    
    this.start=starts
    ############ smoothing parameter loop
    for(spar.i in n.spar:1) { 
#        if(is.infinite(spar[spar.i])){  
#            if(pen.order==1){
#                this.null.fit=scaledTMix.null(tstat,df)
#                parms.new=c(coef(this.null.fit), rep(0,ncol(H)-1))
#                
#                if(tuning.method=='NIC'){
#                    criterion.mean[cv.i, spar.i]=this.null.fit$tuning$mean
#                    criterion.var[cv.i, spar.i]=this.null.fit$tuning$var
#                }else if (tuning.method=='CV'){
#                    pred.ll=logLiktest(parms.new[-1], 1+exp(parms.new[1]))
#                    criterion.mean[cv.i, spar.i]=-mean(pred.ll)
#                    criterion.var[cv.i, spar.i]=var(pred.ll)/length(pred.ll)
#                }else {stop('unknown tuning.method')}
#                enps[cv.i, spar.i]=enps0[cv.i, spar.i]=2
#
#                this.start=parms.new
#                all.parms[,spar.i]=all.parms[,spar.i]+parms.new/cv.fold
#                next
#            }else .NotYetImplemented() ## TODO: Fit global polynomial
#        }
#        if(debugging) cat('############## CV:',cv.i,'\t## spar.i', spar.i,fill=TRUE)
#        spar.Pen.mat=spar[spar.i]*Pen.mat

#        ############## EM iterations
#        if(debugging) time.init=proc.time()[3]
#        em.parms=EMupdate(this.start, nLogLik.pen, optim.method, 
#                          H, tstat, df, dt0, spar.Pen.mat,
#                          em.iter.max, em.beta.iter.max, 
#                          scale.conv, lfdr.conv, NPLL.conv, debugging
#        )
#        if(debugging) cat("\ntime used in EM:", proc.time()[3]-time.init,fill=TRUE)

#        ############## Newtonian iteration
#        if(debugging) time.init=proc.time()[3]
#        if(optim.method%in%c("BFGS","CG","L-BFGS-B","Nelder-Mead", "SANN")){ ## call optim
#            optim.fit=optim(em.parms, nLogLik.pen, deriv.nLogLik.pen, method=optim.method, 
#                        control=list(maxit=newton.iter.max, trace=if(debugging)3 else 0, REPORT=30))
#            parms.new=optim.fit$par
#            if(debugging)cat("optim iter=", optim.fit$iterations,fill=TRUE)
#        }else if (optim.method=='nlminb') {
            #            nlminb.i=1
            #            repeat{
                em.parms=this.start
            nlminb.fit=nlminb(em.parms, nLogLik.pen, deriv.nLogLik.pen, hess.nLogLik.pen, return.dense=TRUE, 
                              control=list(eval.max=newton.iter.max*2, iter.max=newton.iter.max))
            #            if(nlminb.i>=10 || nlminb.fit$convergence==0) break
            #            nlminb.i=nlminb.i+1
            #            em.parms=nlminb.fit$par
            #            if(debugging)cat('another try on nlminb',fill=TRUE)
            #            }
            parms.new=nlminb.fit$par
            if(debugging)cat("nlminb iter=", nlminb.fit$iterations,fill=TRUE)
            if(debugging)cat('max |gradient|=', max(abs(deriv.nLogLik.pen(parms.new))), fill=TRUE)
#        }else if (optim.method=='NR'){
#            NR.fit=NRupdate(nLogLik.pen, em.parms, iter.max=newton.iter.max, with.grad.hess=TRUE,debugging=debugging)
#            parms.new=NR.fit
#        }
#        if(debugging) cat("\ntime used in Newtonian optimization:", proc.time()[3]-time.init,fill=TRUE)



        ############### cross-validation loglik
        pred.ll=logLiktest(parms.new[-1], 1+exp(parms.new[1]))

            ###### enp = tr{ E[Hess]^(-1) Var(Grad) }
#            grad.new=deriv.nLogLik.pen(parms.new,return.K=TRUE)
            hess.new=hess.nLogLik.pen(parms.new)
#            J.inv.K=try(solve(hess.new,attr(grad.new,'K')), silent=TRUE)
#            if(class(J.inv.K)=='try-error'){
#                warning(paste("final Hessian is not positive definite. The smallest eigen value is", tail(eigen(hess.new,TRUE,TRUE)$val,1)))
#                J.inv.K=solve(nearPD(hess.new)$mat, attr(grad.new,'K'))
#            }
#                
#        enp=(sum(diag(J.inv.K)))
        enp=2+n.vars*(pen.order-1)

        ############### CV model selection criteria
        nparm.extra=0       #  1 
        enps[cv.i, spar.i]=enps0[cv.i, spar.i]=enp
        criterion.mean[cv.i, spar.i]=-mean(pred.ll)
        criterion.var[cv.i, spar.i]=var(pred.ll)/length(pred.ll)

        if(debugging){
            cat("enp=",enp,"\tpred.ll=",sum(pred.ll),"\t@spar=",spar[spar.i],fill=TRUE)
        }

        this.start=parms.new
        all.parms[,spar.i]=all.parms[,spar.i]+parms.new/cv.fold
    } ## of smoothing parameter

    ########## logistic correction to the enp
#    if(logistic.correction){
#        enp.logistic=try(logistic.enp(log10(spar), enps[cv.i,], ncol(H.all)+1, 2+n.vars*(pen.order-1)), silent=TRUE)
#        if(class(enp.logistic)=='try-error') enp.logistic=enps[cv.i,]
#        enps[cv.i,]=enp.logistic
#    }
  } ## of cross-validation loop

#        if(tuning.method=='NIC') {
            criterion.mean=criterion.mean +enps/G.all
#        }else if (tuning.method=='CV'){ ## nothing to do here
#            #        }else if (tuning.method=='GCV'){
#            #            criterion.mean[cv.i,spar.i]=sum(pred.deviance)/G/(1-(enp)/G)^2 
#            #            criterion.var[cv.i,spar.i]=var(pred.deviance)/G/(1-(enp)/G)^4
#            #        }else if(tuning.method=='BIC'){
#            #            criterion.mean[cv.i,spar.i]=sum(pred.ll) +log(G)*(enp+nparm.extra)
#            #            criterion.var[cv.i,spar.i]=var(pred.ll)
#            #        }else if(tuning.method=='CAIC'){
#            #            criterion.mean[cv.i,spar.i]=sum(pred.ll) +(log(G)+1)*(enp+nparm.extra)
#            #            criterion.var[cv.i,spar.i]=var(pred.ll)
#            #        }else if(tuning.method=='HQIC'){
#            #            criterion.mean[cv.i,spar.i]=sum(pred.ll) +2*log(log(G))*(enp+nparm.extra)
#            #            criterion.var[cv.i,spar.i]=var(pred.ll)
#        }else{
#            stop("unsupported tuning.method")
#        }


    ############### post-processing
    ## CV is conducted or more than one spar; find good spar and do final call
#    if(!noCV){
#        wt.cv=table(cv.grp)
#        crit.mean.all=crossprod(wt.cv, criterion.mean)/sum(wt.cv)
#        crit.se.all=sqrt(crossprod(wt.cv*wt.cv, criterion.var)/G.all/G.all)
#        s.mode.i=1
#        goodenp.idx=rep(TRUE,n.spar)
#        imin.cv=which.min(crit.mean.all)
#    }else{
        crit.mean.all=drop(criterion.mean)
        crit.se.all=sqrt(drop(criterion.var))
        s.mode.i=1 #if(exists('enp.logistic',inherits=FALSE))which(log10(spar)==attr(enp.logistic,'mode'))[1] else 1
        if(is.na(s.mode.i))s.mode.i=1
        goodenp.idx=rep(FALSE,n.spar)
        goodenp.idx[s.mode.i:n.spar]=TRUE
        imin.cv=s.mode.i-1+which.min(crit.mean.all[goodenp.idx])
#    }


    ############### either final fit on all data or return results
#    if(noCV){ ### directly reture results
        parms.new=all.parms[,imin.cv]
        final.scale=1+exp(parms.new[1])
        final.ll=logLiktest(parms.new[-1], final.scale)
        final.pen=0 #if(is.infinite(spar[imin.cv])) 0 else
            #.5*spar[imin.cv]*drop(parms.new[-1]%*%Pen.mat%*%parms.new[-1])

        fx.final=drop(H.all%*%parms.new[-1])
        pi0.final=1/(1+exp(fx.final))
        dt1.all=dt(tstat.all/final.scale,df)/final.scale
        lfdr.final <- pi0.final*dt0.all
        f.final=lfdr.final+(1-pi0.final)*dt1.all
        lfdr.final=lfdr.final/f.final

        fx.j=matrix(,G,n.vars)
        for(j in 1:n.vars) {j.idx=j==j.all; fx.j[,j]=drop(H.all[,j.idx,drop=FALSE]%*%parms.new[-1][j.idx]); fx.j[,j]=fx.j[,j]-mean(fx.j[,j])}

#        if(is.infinite(spar[imin.cv])){
#            if(!exists('null.fit', inherits=FALSE))null.fit=scaledTMix.null(tstat.all,df)
#            J=K=solve(null.fit$fit$asym.vcov)
            J=K=hess.new
#        }else{
#            J=hess.nLogLik.pen(parms.new)
#            K=attr(deriv.nLogLik.pen(parms.new,TRUE),'K')
#        }
        J.Inv=try(solve(J), silent=TRUE)
        if(class(J.Inv)=='try-error'){
            warning(paste("final Hessian is not positive definite. The smallest eigen value is", tail(eigen(J,TRUE,TRUE)$val,1)))
            J.Inv=symmpart(as(solve(nearPD(J)$mat),'sparseMatrix'))
        }else{
            J.Inv=symmpart(J.Inv)
        }
        cov.parms=J.Inv #symmpart(J.Inv%*%K%*%J.Inv)
        ans=list(lfdr=lfdr.final,
                 model=list(tstat=tstat.all, df=df, x=x, pen.order=pen.order, poly.degree=poly.degree), 
                 scale.fact=list(scale.fact=final.scale, sd.ncp=sqrt(final.scale^2-1),  r=parms.new[1], 
                            t.cross=sqrt(df*(final.scale^(2/(df+1))-1)/(1-final.scale^(-2*df/(df+1))))),
                 pi0=pi0.final,
                 tuning=list(mean=criterion.mean,  var=criterion.var,  grp=cv.grp,  method=tuning.method, final=criterion.mean[imin.cv]),
                 spar=list(all=spar, final=spar[imin.cv], final.idx=imin.cv),
                 enp=list(raw=enps0, logistic=enps, final=enps[,imin.cv], good.idx=goodenp.idx),
                 fit=list(intercept=mean(H.all%*%parms.new[-1]), 
                          covariate.idx=j.all,
                          f.covariate=fx.j,
                          f=fx.final,
                          beta=parms.new[-1],
                          H=H.all,
                          asym.vcov=cov.parms),
                 NPLL=list(NPLL=-sum(final.ll)+final.pen, logLik=final.ll, penalty=final.pen, 
                           saturated.ll=#attr(logLik.saturated,'logLik'),
                                        ifelse(abs(tstat)<sqrt(df*(final.scale^(2/(df+1))-1)/(1-final.scale^(-2*df/(df+1)))), 
                                               dt(tstat,df,log=TRUE), dt(tstat/final.scale,df,log=TRUE)-log(final.scale))
                          )
        )
        class(ans)='hisemit'
        if(plotit)plot(ans)
        return(ans)
#    }
#
#    if(debugging) cat("\n######### FINAL FIT ##############",fill=TRUE)
#
#    final.fit=penLik.EMNewton(tstat.all,x,df, spar=spar[imin.cv], nknots, all.parms[,imin.cv], tuning.method='NIC',
#                cv.fold=1, optim.method, logistic.correction=FALSE,
#                em.iter.max, em.beta.iter.max,newton.iter.max, scale.conv, lfdr.conv,NPLL.conv, debugging, plotit=FALSE)
#    
#    final.fit$tuning=list(mean=criterion.mean,  var=criterion.var,  grp=cv.grp,  method=tuning.method, 
#                    final=final.fit$tuning$final)
#    final.fit$spar=list(all=spar, final=spar[imin.cv], final.idx=imin.cv)
#    final.fit$enp=list(raw=enps0, logistic=enps, final=final.fit$enp$final, good.idx=goodenp.idx)
#    if(plotit)plot(final.fit)
#    return(final.fit)
   
}





if(FALSE){### testing code

    source("scaledTMix.R")
    source("directSum.R")
    source("NRupdate.R")
    source("EMupdate.R")
    source("logistic.enp.R")
    source("EM.Newton.support.R")
    ############################################################################
    ###################            Univariate testing        ###################   
    ############################################################################


    G=20000
    sdncp=1.3
    n1=n2=5
    df=n1+n2-2
    x=1:G
    f=function(x)sin(x*pi/1000)+1
    plot(x,f(x),type='l')
    #Pi.i=rep(.7,G) 
    Pi.i=1/(1+exp(f(x)))
    mean(Pi.i)


    set.seed(97584)
    Z.i=rbinom(G,1,1-Pi.i)
    t0.i=rt(G,df)
    ncp.i=rnorm(G,0,sdncp)
    t1.i=rt(G,df,ncp.i)
    t.i=ifelse(Z.i==0,t0.i,t1.i)
    pvals=2*pt(abs(t.i),df,lower=FALSE)
    require(qvalue); qvalue(pvals)$pi0  
    require(ROCR); performance(prediction(abs(t.i),Z.i),'auc')@y.values[[1]]


    source("scaledTMix.R")
    (stm.null=scaledTMix.null(t.i,df))
    #[1] 0.2461443 1.4399100
    #attr(,"equiv.sd.ncp")
    #[1] 1.036022


    ### true.lfdr
    true.dt1=dt(t.i/sqrt(1+sdncp*sdncp),df)/sqrt(1+sdncp*sdncp)
    true.lfdr <- Pi.i*dt(t.i,df)
    true.ft=true.lfdr+(1-Pi.i)*true.dt1
    true.lfdr=true.lfdr/true.ft
    summary(true.lfdr)
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    #0.004087 0.154900 0.254800 0.290400 0.434900 0.621100 

    source("directSum.R")
    source("NRupdate.R")
    source("EMupdate.R")
    source("min.eigval.R")
    source("EM.Newton.R")
    #sink("tmp.Rout")
    set.seed(23894900)
    tmpcv=penLik.EMNewton(t.i,x,df, tuning.method='CV', spar=10^c(-1:2,seq(2.1,3.9,by=.1),4:8), cv.fold=10, nknots=100, debugging=TRUE)
    #sink(NULL)
    #sink("tmp.Rout")
    tmpnic=penLik.EMNewton(t.i,x,df, tuning.method='NIC', spar=10^c(-1:2,seq(2.1,3.9,by=.1),4:8), cv.fold=1, nknots=100, debugging=TRUE)
    #sink(NULL)

    #png("1d.nic.vs.cv.png")
    plot.cv(tmpcv)
    plot.cv(tmpnic,add=TRUE)
    legend('topright',pch=1:0,lty=0:1,legend=c('CV','NIC'))
    #dev.off()

    tmp=tmpnic  #tmp=tmpcv
    str(tmp)
    #curve(f,1,G, col=2,lwd=3); lines(x,attr(tmp,'H')%*%attr(tmp,'beta'))
    plot(x,attr(tmp,'H')%*%attr(tmp,'beta'),type='l')
    curve(f,1,G, add=TRUE,col=2,lwd=3)
    cor(f(x),drop(attr(tmp,'H')%*%attr(tmp,'beta')))
    cor(true.lfdr, tmp); plot(true.lfdr, tmp); abline(0,1,col=2)
    cor(Pi.i, attr(tmp,'pi0'))


    plot(attr(tmp,'beta'))
    summary(attr(tmp,'pi0'))
    summary(tmp)



    library(ROCR)
    t.pred=prediction(abs(t.i),Z.i)
    ideal.pred=prediction(-true.lfdr,Z.i)
    my.pred=prediction(-tmp[1:G],Z.i)

    performance(t.pred,'auc')@y.values[[1]]
    performance(ideal.pred,'auc')@y.values[[1]]
    performance(my.pred,'auc')@y.values[[1]]

    plot(performance(t.pred,'tpr','fpr'),col=1)
    plot(performance(ideal.pred,'tpr','fpr'),col=4,add=TRUE)
    plot(performance(my.pred,'tpr','fpr'),col=2,add=TRUE)



    ############################################################################
    ###################             Bivariate testing        ###################   
    ############################################################################
    x.null=runif(G)*G
    tmpnic.null=penLik.EMNewton(t.i,cbind(x,x.null),df, tuning.method='NIC', spar=10^c(1:2,seq(2.1,3.9,by=.1),4:8), cv.fold=1, nknots=100, debugging=TRUE)
    tmpcv.null=penLik.EMNewton(t.i,cbind(x,x.null),df, tuning.method='CV', spar=10^c(-1:2,seq(2.1,3.9,by=.1),4:8), cv.fold=10, nknots=100, debugging=TRUE)

    #png("1d.nic.add.null.png")
    plot.cv(tmpnic,ylim=c(1.9225,1.925))
    plot.cv(tmpnic.null,add=TRUE)

    plot.cv(tmpcv,ylim=c(1.9225,1.925),add=TRUE,col=2)
    plot.cv(tmpcv.null,add=TRUE,col=2)
    legend('topright',pch=1:0,lty=0:1,legend=c('NIC: 1D','NIC: 1D+fake', "CV: 1D", "CV: 1D+fake"),col=rep(1:2,each=2))
    #dev.off()


    G=20000
    sdncp=1.3
    n1=n2=5
    df=n1+n2-2
    x=1:G
    set.seed(152356)
    y=sample(x)+runif(G)
    f=function(x,y)sin(x*pi/2000)+1+cos(y*pi/2000)
    #fxy=outer(x,x,f)
    ##library(rgl)
    ##persp3d(x,x,fxy)
    #persp(x,x,fxy)
    #Pi.i=rep(.7,G) 
    Pi.i=1/(1+exp(f(x,y)))
    summary(Pi.i)

    Z.i=rbinom(G,1,1-Pi.i)
    t0.i=rt(G,df)
    ncp.i=rnorm(G,0,sdncp)
    t1.i=rt(G,df,ncp.i)
    t.i=ifelse(Z.i==0,t0.i,t1.i)
    pvals=2*pt(abs(t.i),df,lower=FALSE)
    require(qvalue); qvalue(pvals)$pi0  
    require(ROCR); performance(prediction(abs(t.i),Z.i),'auc')@y.values[[1]]



    source("directSum.R")
    source("NRupdate.R")
    source("EMupdate.R")
    source("min.eigval.R")
    source("EM.Newton.R")
    #sink("tmp.Rout")
    set.seed(23894900)
    tmpcv.2d=penLik.EMNewton(t.i,cbind(x,y),df, tuning.method='CV', spar=10^c(-8:0,seq(1,3,by=.2),4:7), cv.fold=5, nknots=400, debugging=TRUE)
    tmpcv.2d.null=penLik.EMNewton(t.i,(x),df, tuning.method='CV', spar=10^c(-8:0,seq(1,3,by=.2),4:7), cv.fold=5, nknots=400, debugging=TRUE)
    tmpnic.2d=penLik.EMNewton(t.i,cbind(x,y),df, tuning.method='NIC', spar=10^c(-8:0,seq(1,3,by=.2),4:7), cv.fold=10, nknots=400, debugging=TRUE)
    tmpnic.2d.null=penLik.EMNewton(t.i,(x),df, tuning.method='NIC', spar=10^c(-8:1,seq(1,3,by=.1),4,6:7), cv.fold=10, nknots=400, debugging=TRUE)

    #png("2d.diff.cv.nic.png")
    plot.cv(tmpcv.2d,ylim=c(1.928,1.933))
    plot.cv(tmpnic.2d,add=TRUE)
    legend('topright',pch=1:0,lty=0:1,legend=c('CV','NIC'))
    #dev.off()

    #png("2d.nic.drop1d.png")
    plot.cv(tmpnic.2d)
    plot.cv(tmpnic.2d.null,add=TRUE)
    legend('topright',pch=1:0,lty=0:1,legend=c('2D','drop 1D'))
    #dev.off()

    #png("2d.cv.drop1d.png")
    plot.cv(tmpcv.2d,ylim=c(1.928,1.933))
    plot.cv(tmpcv.2d.null,add=TRUE)
    legend('topright',pch=1:0,lty=0:1,legend=c('2D','drop 1D'))
    #dev.off()

    tmp=tmpnic.2d  #tmp=tmpcv.2d    #tmp=tmpnic.2d.null
    str(tmp)
    library(rgl)
    #curve(f,1,G, col=2); lines(x,attr(tmp,'H')%*%attr(tmp,'beta'))
    plot3d(x,y,attr(tmp,'H')%*%attr(tmp,'beta'))
    points3d(x,y,f(x,y),col=2)
    cor(f(x),attr(tmp,'H')%*%attr(tmp,'beta'))
    cor(true.lfdr, tmp); plot(true.lfdr, tmp); abline(0,1,col=2)
    cor(Pi.i, attr(tmp,'pi0'))

    plot.cv(tmp)

    plot(attr(tmp,'beta'))
    summary(attr(tmp,'pi0'))
    summary(tmp)

    source("scaledTMix.R")
    (stm.null=scaledTMix.null(t.i,df))

    ### true.lfdr
    true.dt1=dt(t.i/sqrt(1+sdncp*sdncp),df)/sqrt(1+sdncp*sdncp)
    true.lfdr <- Pi.i*dt(t.i,df)
    true.ft=true.lfdr+(1-Pi.i)*true.dt1
    true.lfdr=true.lfdr/true.ft
    summary(true.lfdr)
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    #0.003288 0.118800 0.263700 0.305200 0.462000 0.811800 
    t.pred=prediction(abs(t.i),Z.i)
    ideal.pred=prediction(-true.lfdr,Z.i)
    my.pred=prediction(-tmp[1:G],Z.i)

    performance(t.pred,'auc')@y.values[[1]]
    performance(ideal.pred,'auc')@y.values[[1]]
    performance(my.pred,'auc')@y.values[[1]]

    plot(performance(t.pred,'tpr','fpr'),col=1)
    plot(performance(ideal.pred,'tpr','fpr'),col=4,add=TRUE)
    plot(performance(my.pred,'tpr','fpr'),col=2,add=TRUE)



}
