
opta=TRUE
aa=opta*poped.db$settings$cfaa*matrix(1,poped.db$design$m,size(poped.db$design$a,2))
aa

optxt=TRUE
axt=optxt*poped.db$settings$cfaxt*matrix(1,poped.db$design$m,max(poped.db$design_space$maxni))
axt

calc_ofv_and_grad(x=c(poped.db$design$xt,poped.db$design$a),
                  optxt=optxt, opta=opta, 
                  model_switch=poped.db$design$model_switch,
                  aa=aa,
                  axt=axt,
                  groupsize=poped.db$design$groupsize,
                  ni=poped.db$design$ni,
                  xtopto=poped.db$design$xt,
                  xopto=poped.db$design$x,
                  aopto=poped.db$design$a,
                  bpop=poped.db$parameters$param.pt.val$bpop,
                  d=poped.db$parameters$param.pt.val$d,
                  sigma=poped.db$parameters$param.pt.val$sigma,
                  docc_full=poped.db$parameters$param.pt.val$docc,
                  poped.db,
                  only_fim=FALSE)

\dontrun{
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)
  
}




