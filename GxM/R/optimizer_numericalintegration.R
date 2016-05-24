optimizer_numericalintegration <-
function(listdata,fMP,zeroset,starting,K,coreNumber) {
   
  npairs = length(listdata);

  ## points and weights in 3d for Gauss-Hermite Quadrature
  PointsWK = ghpoints3(K=K);

  if(coreNumber == 1) {
    OutputK = bobyqa(par=starting$initial,fn=likelihood,
                     lower=starting$lower,upper=starting$upper,control=list(maxfun=3000), 
                     paramnames=starting$paramnames,zeroset=zeroset,listdata=listdata,  
                     fMP=fMP, PointsW=PointsWK);
    estimates = OutputK$par;      
    names(estimates) = starting$paramnames;
    loglike = - OutputK$fval;  
    BIC = -2*loglike + length(estimates)*log(npairs);
    hessfit = fdHess(pars=estimates,fun=likelihood,paramnames=starting$paramnames,zeroset=zeroset,
                     listdata=listdata,fMP=fMP,PointsW=PointsWK);
    hessianmatrix = hessfit$Hessian; colnames(hessianmatrix) = starting$paramnames;
    gradient = hessfit$gradient;
  } ## end of if(coreNumber == 1)


  if(coreNumber > 1) {
    cluster = makeCluster(coreNumber);
    OutputK = bobyqa(par=starting$initial,fn=likelihoodcluster,
                     lower=starting$lower,upper=starting$upper,control=list(maxfun=3000), 
                     paramnames=starting$paramnames,zeroset=zeroset,listdata=listdata,  
                     fMP=fMP,PointsW=PointsWK,cl=cluster);
    estimates = OutputK$par;  names(estimates) = starting$paramnames;
    loglike = - OutputK$fval;  
    BIC = -2*loglike + length(estimates)*log(npairs);
    hessfit = fdHess(pars=estimates,fun=likelihoodcluster,paramnames=starting$paramnames,zeroset=zeroset,
                     listdata=listdata,fMP=fMP,PointsW=PointsWK,cl=cluster);
    hessianmatrix = hessfit$Hessian; colnames(hessianmatrix) = starting$paramnames;
    gradient = hessfit$gradient;
    stopCluster(cluster);
  } ## end of if(coreNumber > 1)

 
  return(list(loglikelihood=loglike, BIC=BIC, par=estimates, hess=hessianmatrix, gradient=gradient));
}
