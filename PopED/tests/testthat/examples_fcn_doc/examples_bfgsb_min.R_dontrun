
\dontrun{  
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=0,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)
  
  f_name <- 'calc_ofv_and_grad' 
  gen_des <- downsizing_general_design(poped.db)
  
  aa <- 0*poped.db$settings$cfaa*matrix(1,poped.db$design$m,size(poped.db$design$a,2))
  axt=1*poped.db$settings$cfaxt*matrix(1,poped.db$design$m,max(poped.db$design_space$maxni))
  
  f_options_1 <- list(gen_des$x,1, 0, gen_des$model_switch,
                    aa=aa,axt=axt,poped.db$design$groupsize,
                    gen_des$ni,
                    gen_des$xt,gen_des$x,gen_des$a,gen_des$bpop[,2,drop=F],
                    getfulld(gen_des$d[,2,drop=F],poped.db$parameters$covd),
                    poped.db$parameters$sigma,
                    getfulld(poped.db$parameters$docc[,2,drop=F],
                             poped.db$parameters$covdocc),poped.db)
  
  options=list('factr'=poped.db$settings$BFGSConvergenceCriteriaMinStep,
               #'factr'=0.01,
               'pgtol'=poped.db$settings$BFGSProjectedGradientTol,
               'ftol'=poped.db$settings$BFGSTolerancef,
               'gtol'=poped.db$settings$BFGSToleranceg,
               'xtol'=poped.db$settings$BFGSTolerancex)
  
  x_k=t(gen_des$xt)
  lb=t(gen_des$minxt)
  ub=t(gen_des$maxxt)
  
  output <- bfgsb_min(f_name,f_options, x_k,lb,ub,options) 
  
}

