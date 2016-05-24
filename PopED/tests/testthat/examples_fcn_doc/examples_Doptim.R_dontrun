
\dontrun{
  
  ##############
  # typically one will use poped_optimize 
  # This then calls Doptim for continuous optimization problems
  ##############
  
  
  # RS+SG+LS optimization of sample times
  # optimization with just a few iterations
  # only to check that things are working
  output <- poped_optimize(poped.db,opt_xt=T,
                           rsit=5,sgit=5,ls_step_size=5)
  
  # RS+SG+LS optimization of sample times 
  # (longer run time than above but more likely to reach a maximum)
  output <- poped_optimize(poped.db,opt_xt=T)
  get_rse(output$fmf,output$poped.db)
  plot_model_prediction(output$poped.db)
  
  
  # Random search (just a few samples here)
  rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                              bUseRandomSearch= 1,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0)
  
  # line search, DOSE and sample time optimization
  ls.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 1,
                              ls_step_size=10)
  
  # Stochastic gradient search, DOSE and sample time optimization
  sg.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1, 
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 1,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0,
                              sgit=20)
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)

  ##############
  # If you really want to you can use Doptim dirtectly
  ##############
  dsl <- downsizing_general_design(poped.db)
  poped.db$settings$optsw[2] <- 1  # sample time optimization
  output <- Doptim(poped.db,dsl$ni, dsl$xt, dsl$model_switch, dsl$x, dsl$a, 
         dsl$bpop, dsl$d, dsl$maxxt, dsl$minxt,dsl$maxa,dsl$mina) 
  
}


