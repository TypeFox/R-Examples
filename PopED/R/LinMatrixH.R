#' Model linearization with respect to epsilon.
#' 
#' The function performs a linearization of the model with respect to the residual variability.
#' Derivative of model w.r.t. eps evaluated at eps=0
#' 
#' @inheritParams mftot
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param xt_ind A vector of the individual/group sample times
#' @param b_ind vector of individual realization of the BSV terms b
#' @param bocc_ind Vector of individual realizations of the BOV terms bocc
#' 
#' @return A matrix of size (samples per individual x number of epsilons) 
#'  
#' @family FIM
#'     
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_LinMatrixH.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixH <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db){
  #----------Model linearization with respect to epsilon.
  #
  # size of return is (samples per individual x number of epsilons) 
  #
  # derivative of model w$r.t. eps eval at e=0
  #
  NumEPS = size(poped.db$parameters$sigma,1)
  if((NumEPS==0)){
    y=0
  } else { 
    returnArgs <- gradf_eps(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,NumEPS,poped.db) 
    y <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
  }
  return(list( y= y,poped.db=poped.db)) 
}

#' Model linearization with respect to epsilon.
#' 
#' The function performs a linearization of the model with respect to the residual variability.
#' Derivative of model w.r.t. eps evaluated at eps=0 and b=b_ind.
#' 
#' @inheritParams mftot
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams LinMatrixH
#' @param num_eps The number of \code{eps()} in the model.
#' 
#' @return A matrix of size (samples per individual x number of epsilons) 
#'  
#' @family FIM
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_gradf_eps.R
#' @export
#' @keywords internal
#' 
gradf_eps <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,num_eps,poped.db){
  #----------Model linearization with respect to epsilon.
  #
  # size of return is (samples per individual x number of epsilons)
  #
  # derivative of model w$r.t. eps eval at e=0 and b=b_ind
  #
  #
  
  dfeps_de0=zeros(size(xt_ind,1),num_eps)
  
  if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==1) ){#No interaction
    fg0=feval(poped.db$model$fg_pointer,x,a,bpop,zeros(size(b_ind)[1],size(b_ind)[2]),zeros(size(bocc_ind)[1],size(bocc_ind)[2]))
    
  } else {
    fg0=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind) #Interaction
  }
  
  e0=zeros(1,num_eps)
  
  #Central approximation
  if((poped.db$settings$hle_switch==1)){
    for(i in 1:num_eps){
      e_plus=e0
      e_minus=e0
      e_plus[i] = e_plus[i]+poped.db$settings$hle
      e_minus[i]= e_minus[i]-poped.db$settings$hle
      returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg0,e_plus,poped.db) 
      ferror_plus <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg0,e_minus,poped.db) 
      ferror_minus <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      dfeps_de0[,i]=(ferror_plus-ferror_minus)/(2*poped.db$settings$hle)
    }
  } else {
    #Complex approximation
    if((poped.db$settings$hle_switch==0)){
      for(i in 1:num_eps){
        e_plus=e0
        e_plus[i] = complex(real=e_plus[i],imaginary=poped.db$settings$hle)
        returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg0,e_plus,poped.db) 
        ferror_plus <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        dfeps_de0[,i]=Im(ferror_plus)/poped.db$settings$hle
      }
    } else {
      if((poped.db$settings$hle_switch==30) ){#Automatic differentiation (INTLab)
        stop("Automatic differentiation not yet implemented in PopED for R")
        #             e_init = gradientinit(e0)
        #              returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg0,e_init,poped.db) 
        # ferror_val <- returnArgs[[1]]
        # poped.db <- returnArgs[[2]]
        #             dfeps_de0 = ferror_val$dx
      } else {
        if((poped.db$settings$hle_switch!=20)){
          stop(sprintf('Unknown derivative option for gradf_eps'))
        }
      }
    }
  }
  return(list( dfeps_de0= dfeps_de0,poped.db=poped.db)) 
}
