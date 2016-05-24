#' Normalize an objective function by the size of the FIM matrix
#' 
#' Compute a normalized OFV based on the size of the FIM matrix.  This value can then be used in 
#' efficiency calculations. This is NOT the OFV used in optimization, see \code{\link{ofv_fim}}. 
#' 
#' @param ofv_f An objective function
#' @param num_parameters The number of parameters to use for normalization
#' @param poped.db a poped database
#' @inheritParams ofv_fim
#' 
#' @return The specified criterion value.
#' 
#' @family FIM
#' 
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_ofv_criterion.R
#' 
#' @export
ofv_criterion <- function(ofv_f,
                          num_parameters,
                          poped.db,
                          ofv_calc_type=poped.db$settings$ofv_calc_type){
  
  #Input: the ofv
  #Return the single value that should be maximized
  
  criterion_value = 0
  
  if((ofv_calc_type==1 || ofv_calc_type==4) ){#D-Optimal Design
    criterion_value = ofv_f^(1/num_parameters)
  }
  
  if((ofv_calc_type==2) ){#A-Optimal Design
    criterion_value=ofv_f/num_parameters
  }
  
  if((ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('Criterion for S-optimal design not implemented yet'))
  }
  
  if((ofv_calc_type==6) ){#Ds-Optimal design
    criterion_value = ofv_f^(1/sum(poped.db$parameters$ds_index))
  }   
  
  return( criterion_value ) 
}
