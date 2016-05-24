
\dontrun{
  
  ##############
  # D-family Optimization
  ##############
  
  # below are a number of ways to optimize the problem
  
  
  # ARS+BFGS+LS optimization of sample times
  # optimization with just a few iterations
  # only to check that things are working
  output <- poped_optim(poped.db,opt_xt=T,
                        control = list(ARS=list(iter=5),
                                       BFGS=list(maxit=5),
                                       LS=list(line_length=5)))
  
  
  # RS+SG+LS optimization of sample times 
  # (longer run time than above but more likely to reach a maximum)
  output <- poped_optim(poped.db,opt_xt=T,parallel = TRUE)
  
  get_rse(output$FIM,output$poped.db)
  plot_model_prediction(output$poped.db)
  
  # optimization with only integer times allowed
  poped.db.2 <- poped.db
  poped.db.2$design_space$xt_space <- matrix(list(seq(1,120)),1,8)
  output_2 <- poped_optim(poped.db.2,opt_xt=T,parallel = TRUE)
  
  get_rse(output_2$FIM,output_2$poped.db)
  plot_model_prediction(output_2$poped.db)
  
  # Examine efficiency of sampling windows
  plot_efficiency_of_windows(output_2$poped.db,xt_windows=0.5)
  plot_efficiency_of_windows(output_2$poped.db,xt_windows=1)
  
  # Random search (just a few samples here)
  rs.output <- poped_optim(poped.db,opt_xt=T,method = "ARS",
                           control = list(ARS=list(iter=5)))
  
  get_rse(rs.output$FIM,rs.output$poped.db)
  
  # line search, DOSE and sample time optimization
  ls.output <- poped_optim(poped.db,opt_xt=T,opt_a=T,method = "LS",
                           control = list(LS=list(line_length=5)))
  
  # Adaptive random search, 
  # DOSE and sample time optimization
  ars.output <- poped_optim(poped.db,opt_xt=T,opt_a=T,method = "ARS",
                           control = list(ARS=list(iter=5)))
  
  # L-BFGS-B gradient search from the stats::optim() function, 
  # DOSE and sample time optimization
  bfgs.output <- poped_optim(poped.db,opt_xt=T,opt_a=T,method = "BFGS",
                            control = list(BFGS=list(maxit=5)))
  
  
  # genetic algorithm from the GA::ga() function, 
  # DOSE and sample time optimization
  ga.output <- poped_optim(poped.db,opt_xt=T,opt_a=T,method = "GA",parallel=T)
  
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
  ars.output <- poped_optim(poped.db,opt_xt=T,opt_a=T,method = "ARS",
                            d_switch=0,use_laplace=TRUE,laplace.fim=TRUE,
                            parallel=T,
                            control = list(ARS=list(iter=5)))
  
}
