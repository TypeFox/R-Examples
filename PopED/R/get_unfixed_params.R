#' Return all the unfixed parameters
#' 
#' all = vector of all unfixed params var derivative is a vector of 1 and 0, 1
#' means derivative of parameter is taken w$r.t. variance otherwise w$r.t. sd If
#' params is supplied then the parameter is taken from this vector instead of
#' poped.db
#' @param poped.db a PopED database.
#' @param params If params is supplied then the parameters are taken from this vector. 
#'   
#' @return A list with the  parameters.  All unfixed parameters are also
#'   returned in the "\code{all} output with the specified order 
#'   (bpop,d,covd,docc,covdocc,sigma,covsigma). \code{var_derivative}  is a
#'   vector of 1's or 0's, 1 means derivative of parameter is taken with respect
#'   to the variance otherwise with respect to standard deviation.
#' @export
#' @keywords internal
get_unfixed_params <- function(poped.db,params=NULL){
  
  if(is.null(params)){
    bpop = poped.db$parameters$bpop[,2,drop=F]
    d = poped.db$parameters$d[,2,drop=F]
    covd = poped.db$parameters$covd
    docc = poped.db$parameters$docc[,2,drop=F]
    covdocc = poped.db$parameters$covdocc
    sigma = diag_matlab(poped.db$parameters$sigma)
    covsigma = zeros(1,length(sigma)*(length(sigma)-1)/2)
    k=1
    for(i in 1:size(poped.db$parameters$sigma,1)){
      for(j in 1:size(poped.db$parameters$sigma,2)){
        if((i<j)){
          covsigma[k] = poped.db$parameters$sigma[i,j]
          k=k+1
        }
      }
    }
  } else {
    nbpop = length(poped.db$parameters$notfixed_bpop)
    nd = length(poped.db$parameters$notfixed_d)
    ncovd = length(poped.db$parameters$notfixed_covd)
    ndocc = length(poped.db$parameters$notfixed_docc)
    ncovdocc = length(poped.db$parameters$notfixed_covdocc)
    nsigma = length(poped.db$parameters$notfixed_sigma)
    ncovsigma = length(poped.db$parameters$notfixed_covsigma)
    
    bpop = params[1:nbpop]
    d=params[(1+nbpop):(nbpop+nd)]
    covd=params[(1+nbpop+nd):(nbpop+nd+ncovd)]
    docc=params[(1+nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc)]
    covdocc=params[(1+nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc)]
    sigma=params[(1+nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma)]
    covsigma=params[(1+nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma)]
  }
  
  bpop = bpop[poped.db$parameters$notfixed_bpop==1]
  d = d[poped.db$parameters$notfixed_d==1]
  covd = covd[poped.db$parameters$notfixed_covd==1]
  docc = docc[poped.db$parameters$notfixed_docc==1]
  covdocc = covdocc[poped.db$parameters$notfixed_covdocc==1]
  sigma = sigma[poped.db$parameters$notfixed_sigma==1]
  covsigma = covsigma[poped.db$parameters$notfixed_covsigma==1]
  
  all = matrix(c(bpop, d, covd, docc, covdocc, sigma, covsigma),ncol=1,byrow=T)
  
  if((poped.db$settings$iFIMCalculationType!=4)){
    var_derivative = matrix(1,size(all))
  } else {
    var_derivative = matrix(c(rep(1,length(bpop)), rep(1,length(d)), rep(1,length(covd)), rep(1,length(docc)), rep(1,length(covdocc)), rep(0,length(sigma)), rep(1,length(covsigma))),ncol=1,byrow=T)
  }
  
  return(list( bpop= bpop,d=d,covd=covd,docc=docc,covdocc=covdocc,sigma=sigma,covsigma=covsigma,all=all,var_derivative =var_derivative )) 
}

