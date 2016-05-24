#MCMC E step for EM collocation
Estep.MCMC = function(times,data,pars,coefs,lik,proc,prior,walk.var,cscale,nstep,in.meth='SplineEst',control.in=NULL)
{
    parhist = matrix(0,nstep,length(pars))                     # Store parameter history
    coefhist = matrix(0,nstep,length(coefs))
    bestcoefhist = matrix(0,nstep,length(coefs))
    acc = rep(0,nstep)

    f = SplineCoefsErr(coefs,times,data,lik,proc,pars,sgn=-1)   # Complete data log posterior

    tcoefs = coefs


    parhist[1,] = pars                              # Initialize
    bestcoefhist[1,]=tcoefs
    coefhist[1,] = tcoefs

    eta = 0        # cscale = 0 defaults to Collocation MCMC

    for(i in 2:nstep){
       print(i)
        if(cscale>0){
          eta = rnorm(length(coefs))*cscale
        }
    tf =  SplineCoefsErr(as.vector(tcoefs)+eta,times,data,lik,proc,pars,sgn=-1)

    print(c(tf,f))
       if( runif(1) < exp( tf-f)){             # Acceptance probability
          coefs = as.vector(tcoefs)+eta
          f = tf
          acc[i] = 1
       }

       parhist[i,] = pars
       coefhist[i,] = coefs
    #   bestcoefhist[i,] = res$coefs
    }

   return(list(parhist=parhist, coefhist=coefhist, acc = acc))
}