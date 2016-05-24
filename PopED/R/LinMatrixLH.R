#' Model linearization with respect to epsilon and eta.
#' 
#' The function performs a linearization of the model with respect to the residual variability and 
#' then the between subject variability.
#' Derivative of model w.r.t. eps then eta, evaluated at eps=0 and b=b_ind.
#' 
#' @inheritParams mftot
#' @inheritParams LinMatrixH
#' @param NumEPS The number of eps() terms in the model.
#' 
#' @return A matrix of size (samples per individual x (number of sigma x number of omega)) 
#'  
#' @family FIM
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_LinMatrixLH.R
#' @export
#' @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

LinMatrixLH <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,NumEPS,poped.db){
  #----------Model linearization with respect to epsilon.
  #
  # size of return is (samples per individual x (number of sigma x number of omega)) 
  #
  # derivative of model w$r.t. sigma then eta, eval at e=0 and eta
  #
  
  
  y = zeros(size(xt_ind,1),poped.db$parameters$NumRanEff*NumEPS)
  if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==1) ){#No interaction
    #return
    return(list( y= y,poped.db=poped.db)) 
  }
  if((poped.db$parameters$NumRanEff==0)){
    #return
    return(list( y= y,poped.db=poped.db)) 
  }
  
  if((poped.db$settings$hle_switch==30) ){#Automatic differentiation (INTLab)
    stop("Automatic differentiation not yet implemented in R version of PopED")
    #     e0=zeros(1,NumEPS)
    #     init_vec = matrix(c(t(b_ind), e0),nrow=1,byrow=T)
    #     deriv_vec = hessianinit(init_vec)
    #      returnArgs <-  new_ferror_file(model_switch,deriv_vec,xt_ind,x,a,bpop,bocc_ind,poped.db) 
    # deriv <- returnArgs[[1]]
    # poped.db <- returnArgs[[2]]
    #     cellDeriv = cell(1,poped.db$parameters$NumRanEff)
    #      for(i in 1:poped.db$parameters$NumRanEff){
    #         tmp=zeros(size(xt_ind,1),NumEPS)
    #         for(j in 1:size(xt_ind,1) ){#for each sample
    #             val = deriv(j)$hx
    #             tmp(j,)=val(poped.db$parameters$NumRanEff+1:poped.db$parameters$NumRanEff+NumEPS,i)
    #         }
    #         cellDeriv[[i]] = tmp
    #      }
    #     
    #      for(i in 1:poped.db$parameters$NumRanEff){
    #         tmp = cellDeriv[[i]]
    #         y[,(i-1)*NumEPS+1:i*NumEPS]=tmp(,1:NumEPS)
    #      }
    #     #return
    #     return(list( y= y,poped.db=poped.db)) 
  }
  
  if((poped.db$settings$hle_switch==20)){
    stop(sprintf('Analytic derivative with interaction is not yet available!'))
  }
  
  for(i in 1:poped.db$parameters$NumRanEff){
    b_ind_plus=b_ind
    b_ind_minus=b_ind
    b_ind_plus[i] = b_ind_plus[i]+poped.db$settings$hle
    b_ind_minus[i]= b_ind_minus[i]-poped.db$settings$hle
    returnArgs <-  LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,poped.db) 
    lin_plus <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    returnArgs <-  LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,poped.db) 
    lin_minus <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    temp = (lin_plus-lin_minus)/(2*poped.db$settings$hle)
    y[,((i-1)*NumEPS+1):(i*NumEPS)]=temp[,1:NumEPS,drop=F]
  }
  return(list( y= y,poped.db=poped.db)) 
}

#Helper function to get the hessian for the AD derivative
new_ferror_file <- function(model_switch,deriv_vec,xt_ind,x,a,bpop,bocc_ind,poped.db){
  fg0=feval(poped.db$model$fg_pointer,x,a,bpop,deriv_vec(1:poped.db$parameters$NumRanEff),bocc_ind) #Interaction
  returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg0,deriv_vec(poped.db$parameters$NumRanEff+1:length(deriv_vec)),poped.db) 
  f_error <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  return(list( f_error= f_error,poped.db=poped.db)) 
}
