
\dontrun{
  
  ##############
  # D-family Optimization
  ##############
  
  # below are a number of ways to optimize the problem
  
  
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
  
  # MFEA optimization with only integer times allowed
  mfea.output <- poped_optimize(poped.db,opt_xt=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(mfea.output$fmf,mfea.output$poped.db)
  plot_model_prediction(mfea.output$poped.db)
  
  # Examine efficiency of sampling windows
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5)
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1)
  
  # Random search (just a few samples here)
  rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                              bUseRandomSearch= 1,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0)
  get_rse(rs.output$fmf,rs.output$poped.db)
  
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
  # E-family Optimization
  ##############
  
  # Adding 10% log-normal Uncertainty to fixed effects (not Favail)
  bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
  bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                           bpop_vals,
                           ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
  bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
  bpop_vals_ed_ln
  
  ## -- Define initial design  and design space
  poped.db <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                    fg_file="sfg",
                                    fError_file="feps.add.prop",
                                    bpop=bpop_vals_ed_ln, 
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6), 
                                    sigma=c(0.01,0.25),
                                    groupsize=32,
                                    xt=c( 0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=120,
                                    a=70,
                                    mina=0,
                                    maxa=100)
  
  # ED optimization using Random search (just a few samples here)
  output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=10,d_switch=0)
  get_rse(output$fmf,output$poped.db)
  
  # ED with laplace approximation, 
  # optimization using Random search (just a few samples here)
  output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=10,
                           d_switch=0,use_laplace=TRUE,laplace.fim=TRUE)
  get_rse(output$fmf,output$poped.db)
  
  
}
