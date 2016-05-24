optimizer_closedform <-
function(listdata,fMP,zeroset,starting,coreNumber) {

  npairs = length(listdata);
  
  if(coreNumber == 1) {
    Output = bobyqa(par=starting$initial,fn=likelihood_closedform,
                    lower=starting$lower,upper=starting$upper,control=list(maxfun=3000), 
                    paramnames=starting$paramnames,zeroset=zeroset,listdata=listdata,fMP=fMP);
    estimates = Output$par;  names(estimates) = starting$paramnames;
    loglike = - Output$fval;  
    BIC = -2*loglike + length(estimates)*log(npairs);
    hessfit = fdHess(pars=estimates,fun=likelihood_closedform,paramnames=starting$paramnames,zeroset=zeroset,listdata=listdata,fMP=fMP);
    hessianmatrix = hessfit$Hessian; colnames(hessianmatrix) = starting$paramnames;
    gradient = hessfit$gradient;
  } ## end of if(coreNumber==1)

  if(coreNumber > 1) {
    cluster = makeCluster(coreNumber);
    Output = bobyqa(par=starting$initial,fn=likelihoodcluster_closedform,
                    lower=starting$lower,upper=starting$upper,control=list(maxfun=3000), 
                    paramnames=starting$paramnames,zeroset=zeroset,listdata=listdata,fMP=fMP,cl=cluster);
    estimates = Output$par;   names(estimates) = starting$paramnames;
    loglike = - Output$fval;  
    BIC = -2*loglike + length(estimates)*log(npairs);
    hessfit = fdHess(pars=estimates,fun=likelihoodcluster_closedform,paramnames=starting$paramnames,zeroset=zeroset,listdata=listdata,fMP=fMP,cl=cluster);
    hessianmatrix = hessfit$Hessian; colnames(hessianmatrix) = starting$paramnames;
    gradient = hessfit$gradient;
    stopCluster(cluster);
  } ## end of if(coreNumber>1)

  return(list(loglikelihood=loglike, BIC=BIC, par=estimates, hess=hessianmatrix, gradient=gradient));
}
