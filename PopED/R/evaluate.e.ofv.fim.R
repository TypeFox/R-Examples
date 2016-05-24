#' Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).
#' 
#' Compute the expectation of the FIM and OFV(FIM) given the model, parameters, distributions of parameter uncertainty, design and methods defined in the 
#' PopED database. Some of the arguments coming from the PopED database can be overwritten;  
#' by default these arguments are \code{NULL} in the 
#' function, if they are supplied then they are used instead of the arguments from the PopED database.
#' 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param use_laplace Should the Laplace method be used in calculating the expectation of the OFV?  
#' @param laplace.fim Should an E(FIM) be calculated when computing the Laplace approximated E(OFV).  Typically
#' the FIM does not need to be computed and, if desired,  this calculation
#' is done usng the standard MC integration technique, so can be slow. 
#' 
#' @return A list containing the E(FIM) and E(OFV(FIM)) and the a poped.db updated according  to the function arguments.
#' 
#' @family FIM
#' @family E-family
#' @family evaluate_FIM
#'  
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_ed.R
#' @example tests/testthat/examples_fcn_doc/examples_evaluate.e.ofv.fim.R
#' @export


evaluate.e.ofv.fim <- function(poped.db,
                               fim.calc.type=NULL,
                               bpop=poped.db$parameters$bpop,
                               d=poped.db$parameters$d,
                               covd=poped.db$parameters$covd,
                               docc=poped.db$parameters$docc,
                               sigma=poped.db$parameters$sigma,
                               model_switch=NULL,
                               ni=NULL,
                               xt=NULL,
                               x=NULL,
                               a=NULL,
                               groupsize=poped.db$design$groupsize,
                               deriv.type = NULL,
                               bLHS=poped.db$settings$bLHS,
                               ofv_calc_type = poped.db$settings$ofv_calc_type,
                               ED_samp_size = poped.db$settings$ED_samp_size,
                               use_laplace=poped.db$settings$iEDCalculationType, 
                               laplace.fim=FALSE, 
                               ...){

  ## update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
      
    }
  }
  
  downsize.list <- downsizing_general_design(poped.db)
  if(is.null(ni)) ni <- downsize.list$ni
  if(is.null(xt)) xt <- downsize.list$xt
  if(is.null(model_switch)) model_switch <- downsize.list$model_switch
  if(is.null(x)) x <- downsize.list$x
  if(is.null(a)) a <- downsize.list$a    
  
  if(is.null(groupsize)) groupsize <- poped.db$design$groupsize
  
  if(!is.null(fim.calc.type)) poped.db$settings$iFIMCalculationType=fim.calc.type
  
  if(!is.null(deriv.type)){ 
    poped.db$settings$m1_switch=deriv.type
    poped.db$settings$m2_switch=deriv.type
    poped.db$settings$hle_switch=deriv.type
    poped.db$settings$gradff_switch=deriv.type
    poped.db$settings$gradfg_switch=deriv.type
  }
  
  E_fim <- NULL
  E_ofv <- NULL
  
  if(!use_laplace){
    output <- ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped.db)
    E_fim <- output$ED_fim
    E_ofv <- output$ED_ofv
    poped.db=output$poped.db
  } else { 
    #stop("Laplce method not yet implemented in R version of PopED")
    E_ofv  <- ed_laplace_ofv(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped.db,...)[["f"]]   
    if(laplace.fim) {
      E_fim <- ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped.db)[["ED_fim"]]
    }
  }    
  return(list(E_ofv=E_ofv,E_fim= E_fim, poped.db=poped.db))
}
