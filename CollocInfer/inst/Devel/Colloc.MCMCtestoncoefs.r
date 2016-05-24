Colloc.MCMC = function(times,data,pars,coefs,lik,proc,prior,walk.var,cscale,nstep,in.meth='SplineEst',control.in=NULL)
{
    parhist = matrix(0,nstep,length(pars))                     # Store parameter history
    coefhist = matrix(0,nstep,length(coefs))
    bestcoefhist = matrix(0,nstep,length(coefs))
    acc = rep(0,nstep)
    
    ei = eigen(walk.var)                                       # Random walk covariance
    S = ei$vectors%*%diag(sqrt(ei$values))

     pr = prior(pars)
    
    f = SplineCoefsErr(coefs,times,data,lik,proc,pars,sgn=-1)  + pr  # Complete data log posterior
#    H = SplineCoefsDC2(coefs,times,data,lik,proc,pars,sgn=-1)
#    G = SplineCoefsDCDP(coefs,times,data,lik,proc,pars,sgn=-1)
    
#    dcdp = ginv(H)%*%G
 
    tcoefs = coefs
  #  
#     ei = eigen(tH)
#     eta = rnorm(length(coefs))
#     logq = sum(eta^2)/2   
#
#     coefs = coefs + ei$vectors%*%diag(1/sqrt(ei$values*(ei$values>0)))%*%eta            
# 

    parhist[1,] = pars                              # Initialize
    bestcoefhist[1,]=tcoefs
    coefhist[1,] = tcoefs

    eta = 0        # cscale = 0 defaults to Collocation MCMC
    logq = 0
    for(i in 2:nstep){
       print(c(i,pars))

       delta = S %*% rnorm(length(pars))   # Random walk step
       
       # Find the nearest parameters in the history

       dists = apply( (parhist[1:(i-1),] - matrix(pars+delta,i-1,length(pars),byrow=TRUE))^2,1,sum)
       whichrow = sort( dists,decreasing=FALSE,index.return=TRUE)$ix[1]
       tcoefs = bestcoefhist[whichrow,]
       
       # Minimization
         
       res = inneropt(times=times,data=data,coefs=matrix(as.vector(tcoefs), nrow=ncol(lik$bvals),length(coefs)/ncol(lik$bvals)),lik=lik,proc=proc,
            pars=pars+delta,in.meth=in.meth,control.in=control.in)
    
       tcoefs = res$coefs                                                  # Now a Gaussian error about the minimum
       tH = SplineCoefsDC2(tcoefs,times,data,lik,proc,pars+delta,sgn=-1)
         
       tpr = prior(pars+delta)
       if(cscale>0){
         ei = eigen(tH)
         eta = rnorm(length(coefs))
  #       logq = -sum(eta^2)/2 - sum(log(ei$values[ei$values>0])/2)
          logq = sum(dnorm(eta,log=TRUE))
         
         eta = cscale*ei$vectors%*%diag(sqrt(ei$values*(ei$values>0)))%*%eta
       }
       
       tf =  SplineCoefsErr(as.vector(tcoefs)+eta,times,data,lik,proc,pars+delta,sgn=-1) + tpr - logq
       
print(c(tf,f,logq,pars+delta))
       if( runif(1) < exp( tf-f)){             # Acceptance probability
          coefs = as.vector(tcoefs)+eta
          pars = pars + delta
          pr = tpr
          
          f = tf
#          H = SplineCoefsDC2(tcoefs,times,data,lik,proc,pars,sgn=-1)
#          G = SplineCoefsDCDP(tcoefs,times,data,lik,proc,pars,sgn=-1)
#          dcdp = ginv(H)%*%G

          acc[i] = 1
       }
       
       parhist[i,] = pars
       coefhist[i,] = coefs
       bestcoefhist[i,] = res$coefs
    }

   return(list(parhist=parhist, coefhist=coefhist, acc = acc))
}