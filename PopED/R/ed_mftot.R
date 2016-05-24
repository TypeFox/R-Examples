#' Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).
#' 
#' Compute the expectation of the FIM given the model, parameters, distributions of parameter uncertainty, design and methods defined in the 
#' PopED database. 
#' 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param xtoptn The xtoptn value
#' @param xoptn The xoptn
#' @param aoptn The aoptn value
#' 
#' @return A list containing the E(FIM) and E(OFV(FIM)) and the a poped.db.
#' 
#' @family FIM
#' @family E-family
#'  
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_ed.R
#' @example tests/testthat/examples_fcn_doc/examples_ed_mftot.R
#' @export
#' @keywords internal
#' 
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

ed_mftot <- function(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped.db){
  #+++++++++++++++++++++++ ED OFV(MF) VALUE
  s=0
  s1=0
  
  fim_list=cell(1,poped.db$settings$ED_samp_size)
  d_gen_list=cell(1,poped.db$settings$ED_samp_size)
  docc_gen_list=cell(1,poped.db$settings$ED_samp_size)
  
  
  bpop_gen  <-  pargen(bpopdescr,poped.db$model$user_distribution_pointer,
                    poped.db$settings$ED_samp_size,poped.db$settings$bLHS,zeros(1,0),poped.db)
  
  for(ct in 1:poped.db$settings$ED_samp_size){
    d_gen = getfulld(pargen(ddescr,poped.db$model$user_distribution_pointer,1,poped.db$settings$bLHS,ct,poped.db),covd)
    docc_gen = getfulld(pargen(docc,poped.db$model$user_distribution_pointer,1,poped.db$settings$bLHS,ct,poped.db),poped.db$parameters$covdocc)
    returnArgs <- mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop_gen[ct,],d_gen,sigma,docc_gen,poped.db) 
    mftmp <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    s=s+ofv_fim(mftmp,poped.db)
    s1=s1+mftmp
    fim_list[[ct]]=mftmp
    d_gen_list[[ct]]=d_gen
    docc_gen_list[[ct]]=docc_gen
  }
  if((!isempty(poped.db$settings$ed_penalty_pointer))){
    returnArgs <- feval(poped.db$settings$ed_penalty_pointer,fim_list,bpop_gen,d_gen_list,docc_gen_list,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped.db) 
    ED_fim <- returnArgs[[1]]
    ED_ofv <- returnArgs[[2]]
    poped.db <- returnArgs[[3]]
  } else {
    ED_ofv=s/poped.db$settings$ED_samp_size
    ED_fim=s1/poped.db$settings$ED_samp_size
  }
  return(list( ED_fim= ED_fim,ED_ofv=ED_ofv,poped.db=poped.db)) 
}

