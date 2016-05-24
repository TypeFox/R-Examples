#' Evaluate a criterion of the Fisher Information Matrix (FIM)
#' 
#' Compute a criterion of the FIM given the model, parameters, design and methods defined in the 
#' PopED database. 
#' 
#' @param fmf The FIM
#' @param poped.db A poped database
#' @param ofv_calc_type  OFV calculation type for FIM
#' \itemize{ 
#' \item 1 = "D-optimality". Determinant of the FIM: det(FIM)
#' \item 2 = "A-optimality".  Inverse of the sum of the expected parameter variances: 
#' 1/trace_matrix(inv(FIM)) 
#' \item 4 = "lnD-optimality".  Natural logarithm of the determinant of the FIM: log(det(FIM)) 
#' \item 6 = "Ds-optimality". Ratio of the Determinant of the FIM and the Determinant of the uninteresting
#' rows and columns of the FIM: det(FIM)/det(FIM_u)
#' \item 7 = Inverse of the sum of the expected parameter RSE: 1/sum(get_rse(FIM,poped.db,use_percent=FALSE))
#' }
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' 
#' @return The specified criterion value.
#' 
#' @family FIM
#' @family evaluate_FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_ofv_fim.R
#' @export
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

ofv_fim <- function(fmf,poped.db,
                    ofv_calc_type = poped.db$settings$ofv_calc_type,
                    ds_index=poped.db$parameters$ds_index,...){
  
  #Input: the FIM
  #Return the single value that should be maximized
  
  ## create ds_index vector if not already done
  if(!is.matrix(ds_index)) ds_index <- matrix(ds_index,1,length(ds_index))
  
  ofv_value = 0
  
  if((!isempty(poped.db$settings$prior_fim) && all(size(poped.db$settings$prior_fim)==size(fmf)))){
    fmf = fmf + poped.db$settings$prior_fim
  }
  
  if((ofv_calc_type==1) ){#determinant of FIM
    if((isempty(fmf))){
      ofv_value = 0
    } else {
      ofv_value = det(fmf)
    }
    return
  }
  
  if((ofv_calc_type==2) ){#trace of the inverse of FIM
    imf = inv(fmf)
    ofv_value = trace_matrix(imf)
    ofv_value = 1/ofv_value #Make it a max-problem
    return
  }
  
  if((ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('S-optimal design not implemented yet!'))
  }
  
  if((ofv_calc_type==4) ){#log determinant of FIM
    #ofv_value = sum(log(svd(fmf)))
    det_fim <- det(fmf)
    if(det_fim<0) det_fim <- 0
    ofv_value = log(det_fim)
  }
  
  if((ofv_calc_type==5) ){# C-optimal design
    stop(sprintf('C-optimal design is not implemented yet!'))
  }
  if((ofv_calc_type==6) ){#Ds-optimal design
    tmp = fmf
    tmp <- tmp[c(col(ds_index)[ds_index==1]),,drop=F]
    tmp <- tmp[,c(col(ds_index)[ds_index==1]),drop=F]
    ofv_value = det(fmf)/det(tmp)
  }
  if((ofv_calc_type==7) ){#sum of CV
    if((sum(sum(fmf))!=0 && !isnan(sum(sum(fmf))))){
      imf = inv(fmf)
      returnArgs <-  get_cv(diag_matlab(imf),poped.db$parameters$bpop,poped.db$parameters$d,poped.db$parameters$docc,poped.db$parameters$sigma,poped.db) 
      params <- returnArgs[[1]]
      params_cvs <- returnArgs[[2]]
      if((isnan(sum(diag_matlab(imf))))){
        ofv_value = 0
      } else {
        ofv_value=1/sum(abs(params_cvs))
      }
    } else {
      ofv_value = 0
    }
    return
  }
  return( ofv_value ) 
}
