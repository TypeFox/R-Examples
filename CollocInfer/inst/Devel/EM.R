Qfn = function(pars,times,data,coefs,lik,proc,oldpars,tfn=NULL,active=1:length(pars))   # EM algorithm Q function
{
   
   cpars = oldpars
   cpars[active] = pars

    if( !is.null(tfn)){
    cpars = tfn(cpars)
    oldpars = tfn(oldpars)
    }

    f = SplineCoefsErr(coefs,times,data,lik,proc,cpars)

    H1 = SplineCoefsDC2(coefs,times,data,lik,proc,cpars)    
    V = SplineCoefsDC2(coefs,times,data,lik,proc,oldpars)

   F = f + 0.5*sum(diag(H1%*%ginv(V)))

    return( F )
}



EM = function(pars,times,data,coefs,lik,proc,niter=10,in.meth,control.in,tfn=NULL,active=1:length(pars))  # EM algorithm
{
    tpars = pars[active]
 
    for(i in 1:niter){

      pars2 = pars
      if(!is.null(tfn)){ 
        pars2 = tfn(pars) 
        print(pars2)
      } 

#       Cres = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
#        control=control.in,times=times,data=data,lik=lik,proc=proc,pars=pars2)

#        coefs = Cres$par
        Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
        coefs = Ires$coefs
        Cres = Ires$res
    
        Pres = nlminb(tpars,Qfn,gr=NULL,hessian=T,control=control.in,
            times=times,data=data,coefs=coefs,lik=lik,proc=proc,oldpars=pars,tfn=tfn,active=active)

        tpars = Pres$par

        pars[active] = tpars

        print(c(i,Pres$obj,pars))
    }

    return( list(pars,coefs,Pres=Pres,Cres=Cres) )
}
