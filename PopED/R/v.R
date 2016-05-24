## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

v <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db){
  # number of samples X number of samples (per individual)
  
  bUseAutoCorrelation = !isempty(poped.db$model$auto_pointer)
  
  bUseFullSigmaCorrelation = FALSE
  
  if((poped.db$settings$m2_switch[1]==0 || poped.db$settings$m2_switch[1]==1 || poped.db$settings$m2_switch[1]==30)){
    returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db) 
    l <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    returnArgs <- LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db) 
    h <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    ret = zeros(0,1)
    
    if((!isempty(sigma) && bUseFullSigmaCorrelation) ){#Update sigma to be fully correlated
      for(i in 1:size(sigma,1)){
        for(j in 1:size(sigma,2)){
          if((i!=j)){
            sigma[i,j]=sqrt(sigma[i,i]*sigma[j,j])
          }
        }
      }
    }
    
    #Add all IIV
    if((!isempty(d))){
      ret = l%*%d%*%t(l)
    } else {
      ret = zeros(length(xt_ind), length(xt_ind))
    }
    
    if((poped.db$settings$bUseSecondOrder)){
      var_eta = zeros(1,length(xt_ind))
      for(o in 1:length(xt_ind)){
        hessian_eta = hessian_eta_complex(model_switch,xt_ind[o],x,a,bpop,b_ind,bocc_ind,poped.db)
        var_eta[o] = 1/4*trace_matrix(hessian_eta*d*(2*hessian_eta)*d)
      }
      ret = ret+diag_matlab(var_eta)
    }
    
    locc = cell(1,poped.db$parameters$NumOcc)
    
    #Add all occasion variability
    for(i in 1:poped.db$parameters$NumOcc){
      if(poped.db$parameters$NumOcc==0) next
      returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,i,poped.db) 
      locc_tmp <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      if((isempty(ret))){
        ret = locc_tmp%*%docc%*%t(locc_tmp)
      } else {
        ret = ret+locc_tmp%*%docc%*%t(locc_tmp)
      }
      locc[[i]] = locc_tmp
    }
    
    if((!isempty(sigma)) ){#If we have any residual variance
      interact = zeros(0,1) 
      if((!isempty(d)) ){#Calculate the interaction terms
        returnArgs <- LinMatrixLH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,size(sigma,1),poped.db) 
        lh <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        interact=lh%*%kron_tmp(d,sigma)%*%t(lh)
        ret = ret + diag_matlab(diag_matlab(interact))
      }
      
      if((bUseAutoCorrelation) ){#Add autocorrelation
        autocorr = feval(poped.db$model$auto_pointer,h,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,l,locc,interact,poped.db)
        if((poped.db$settings$m2_switch[1]==30)){
          stop(sprintf('User definied variance structure could not be used with automatic differentiation'))
        }
        if((isempty(ret))){
          ret = autocorr
        } else {
          ret = ret+ autocorr
        }
      } else { #Add linearized residual model
        if((isempty(ret))){
          ret = diag_matlab(diag_matlab(h%*%sigma%*%t(h)))
        } else {
          ret = ret + diag_matlab(diag_matlab(h%*%sigma%*%t(h)))
        }
      }
    }
  } else {
    if((poped.db$settings$m2_switch[1]==20)){
      stop("Analytic variance not defined")
      #ret = analytic_variance(xt_ind,x,a,bpop,d)
    } else {
      stop(sprintf('Unknown derivative option for variance'))
    }
  }
  return(list( ret= ret,poped.db =poped.db )) 
}

## %This function fixes a bug in FreeMat 4.0
## function ret=diag(a)
##     if (~isempty(a) && size(a,1)==1 && size(a,2)==1)
##         ret=builtin('diag',[a]);
##     else
##         ret=builtin('diag',a);
##     end
## end

