## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error to 
##   avoid sample times at very low concentrations (time 0 or very late samoples).
library(PopED)


# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 

iNumSimulations <- ifelse(fast,5,100)
EAStepSize <- ifelse(fast,40,1)
rsit <- ifelse(fast,3,300)
sgit <- ifelse(fast,3,150)
ls_step_size <- ifelse(fast,3,50)
iter_max <- ifelse(fast,1,10)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list(y=y,poped.db=poped.db))
  })
}

feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional + additive
  y <- ff(model_switch,xt,parameters,poped.db)[[1]] 
  y = y*(1+epsi[,1]) + epsi[,2]  
  return(list(y=y,poped.db=poped.db)) 
}


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff",
                                  fg_file="sfg",
                                  fError_file="feps",
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
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

##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability 
plot_model_prediction(poped.db,IPRED=T,DV=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

##############
# Optimization
##############

# below are a number of ways to optimize the problem

# RS+SG+LS optimization of sample times
output <- poped_optimize(poped.db,opt_xt=T,
                         rsit=rsit,sgit=sgit,ls_step_size=ls_step_size,
                         iter_max=iter_max)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db)


# MFEA optimization with only integer times (or steps of 40 units if fast==TRUE) allowed
mfea.output <- poped_optimize(poped.db,opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=EAStepSize)
get_rse(mfea.output$fmf,mfea.output$poped.db)
plot_model_prediction(mfea.output$poped.db)

# Examine efficiency of sampling windows
plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5,iNumSimulations = iNumSimulations)


# Random search optimization of DOSE and sampling times
rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=rsit)
rs.output$xt
names(rs.output)
result.db <- rs.output$poped.db
plot_model_prediction(result.db)

# RS using the full FIM
rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=rsit,fim.calc.type=0) 

# RS within poped_optimize 
rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                            rsit=rsit,
                            bUseRandomSearch= 1,
                            bUseStochasticGradient = 0,
                            bUseBFGSMinimizer = 0,
                            bUseLineSearch = 0)
names(rs.output)
get_rse(rs.output$fmf,rs.output$poped.db)

# line search, DOSE and sample time optimization
ls.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                            bUseRandomSearch= 0,
                            bUseStochasticGradient = 0,
                            bUseBFGSMinimizer = 0,
                            bUseLineSearch = 1,
                            ls_step_size=ls_step_size,
                            iter_max=iter_max)


# Stochastic gradient search, DOSE and sample time optimization
sg.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1, 
                            bUseRandomSearch= 0,bUseStochasticGradient = 1,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            sgit=sgit)

# BFGS search, DOSE and sample time optimization
if(fast) poped.db$settings$BFGSConvergenceCriteriaMinStep <- 0.01
bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                              bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 1,bUseLineSearch = 0)

