smooth.skeleton <- function(pomp,data,times,params,coefs,basisvals,lambda,in.meth,control.in,eps=1e-6)
{

    profile.obj = setup.profile(pomp,coefs,basisvals,lambda,params,eps)

    lik = profile.obj$lik
    proc= profile.obj$proc


    if(in.meth=="house"){
        res = SplineEst.NewtRaph(coefs,times,data,lik,proc,params,control.in)
        ncoefs = matrix(res$coefs,ncol(lik$bvals),length(res$coefs)/ncol(lik$bvals))
    }
    else if(in.meth=="BFGS"){
        res = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=TRUE,
            control=control.in,times=times,data=data,lik=lik,proc=proc,pars=params,method="BFGS")

        ncoefs = matrix(res$par,ncol(lik$bvals),length(res$par)/ncol(lik$bvals))
    }    
    else if(in.meth=="nlminb"){
        res = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
            control=control.in,times=times,data=data,lik=lik,proc=proc,pars=params)

        ncoefs = matrix(res$par,ncol(lik$bvals),length(res$par)/ncol(lik$bvals))
    }
    else if(in.meth=="maxNR"){
        res = maxNR(SplineCoefsErr,coefs,times=times,data=data,lik=lik,proc=proc,pars=params,sgn=-1,
            grad=SplineCoefsDC,hess=SplineCoefsDC2,print.level=control.in$print.level,
            iterlim = control.in$iterlim)

        ncoefs = matrix(res$estimate,ncol(lik$bvals),length(res$estimate)/ncol(lik$bvals))
    }

    attr(ncoefs,'lik') = lik
    attr(ncoefs,'proc') = proc

    return( ncoefs )

}


profile.skeleton <- function(pomp,data,times,params,coefs,basisvals,lambda,
                                in.meth,out.meth,control.in,control.out,eps=1e-6,active=NULL,covariance=FALSE)
{
    if(is.null(active)){ active = 1:length(params) }

    profile.obj = setup.profile(pomp,coefs,basisvals,lambda,params,eps)

    lik = profile.obj$lik
    proc= profile.obj$proc

#    ProfileEnv = new.env()
#    assign('optcoefs',coefs,3,ProfileEnv)
#    assign('curcoefs',coefs,3,ProfileEnv)

    write.table(coefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
    write.table(coefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)
    write.table(0,file='counter.tmp',col.names=FALSE,row.names=FALSE)
    
    aparams = params[active]
    aparamnames = names(aparams)

    if(out.meth == "BFGS"){
        res = optim(aparams,ProfileErr,pars=params,times=t,data=data,coefs=coefs,lik=lik,proc=proc,hessian=T,
           in.meth=in.meth,control.in=control.in,active=active,control=control.out,gr=ProfileDP.Active,method="BFGS")
        aparams = res$par
    }
    else if(out.meth == "nlminb"){    
        res = nlminb(aparams,ProfileErr,pars=params,times=t,data=data,coefs=coefs,lik=lik,proc=proc,
           in.meth=in.meth,control.in=control.in,active=active,control=control.out,gr=ProfileDP.Active)
        aparams = res$par
    }
    else if(out.meth == "maxNR"){
        require('maxLik')
        res = maxNR(ProfileErr,start=aparams,pars=params,times=t,data=data3,coefs=coefs,lik=lik3,proc=proc,
            in.meth=in.meth,control.in=control.in,active=active,sgn=-1,grad=ProfileDP.Active,
            print.level=control.out$print.level,iterlim=control.out$iterlim,active.par = control.out$activepar)
        aparams = res$est
    }

    names(aparams) = aparamnames

    params[active] = aparams

    #ncoefs = get('optcoefs',envir=ProfileEnv)

    ncoefs = as.matrix(read.table(file='curcoefs.tmp'))
    
    Covar = NULL
    
    if(covariance){
        Covar = profile.covariance(params,active,t,data,ncoefs,lik,proc)
    }


    file.remove('counter.tmp')
    file.remove('curcoefs.tmp')
    file.remove('optcoefs.tmp')

    return( list( params=params, Cov=Covar, res=res, lik=lik, proc=proc ) )
}



setup.profile = function(pomp,coefs,basisvals,lambda,params,eps=1e-6)
{
    lik = make.findif.dmeasure()
    lik$more = list(pomp=pomp,eps=eps)    

    lik$bvals = basisvals$basis.obs
    
    lik$more$names = colnames(coefs)
    lik$more$parnames = names(params)

    proc = make.SSEproc()
    proc$bvals = list(bvals=basisvals$bvals,dbvals=basisvals$dbvals)
    proc$more = make.findif.skeleton()


    if(is.basis(basisvals)){
      if(is.null(quadrature) & is.null(times)){
        error('if basisvals is is a basis object, you must specify the observation times')
      }

      lik$bvals = eval.basis(times,basisvals)

      qpts = basisvals$params[-length(basisvals$params)] + diff(basisvals$params)

      proc$bvals = list()
      proc$bvals$bvals = eval.basis(qpts,basisvals)
      proc$bvals$dbvals = eval.basis(qpts,basisvals,1)

      proc$more$weights = rep(1,length(qpts))
      
    }
    else{

      lik$bvals = basisvals$bvals.obs

      proc$bvals =  list(bvals=basisvals$bvals,dbvals=basisvals$dbvals)
      proc$more$qpts = basisvals$qpts
      proc$more$weights = basisvals$qwts
      if(is.null(proc$more$weights)){
          proc$more$weights = array(1,c(length(proc$more$qpts),dim(coefs)[2]))    
      }
    }
 
    if(length(lambda)==1){ lambda = rep(lambda,dim(coefs,2)) }
    
    
    proc$more$weights = proc$more$weights%*%diag(lambda)

    proc$more$names = colnames(coefs)
    proc$more$parnames = names(params)
    
    proc$more$more = list(pomp=pomp,eps=eps)    

    return( list(lik=lik,proc=proc) )

}
