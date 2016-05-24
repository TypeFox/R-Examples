## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error to 
##   avoid sample times at very low concentrations (time 0 or very late samoples).

## discrete optimization

# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 

EAStepSize <- ifelse(fast,40,1)
rsit <- ifelse(fast,3,300)
ls_step_size <- ifelse(fast,3,50)

########################
## method 1: use the MFEA method with a unit stepsize
########################
library(PopED)

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

# MFEA optimization with only integers (or multiples of 40 if fast=TRUE) 
# in xt and dose allowed (or original design)
# faster optimization than RS+SG+LS in this case
mfea.output <- poped_optimize(poped.db,opt_xt=T,opt_a=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=EAStepSize)
get_rse(mfea.output$fmf,mfea.output$poped.db)
plot_model_prediction(mfea.output$poped.db)


########################
## method 2: use discrete design variables
########################

sfg.discrete <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( CL=bpop[1]*exp(b[1]),
                V=bpop[2]*exp(b[2]),
                KA=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=x[1])
  return( parameters ) 
}

x.space <- cell(1,1)
x.space[1,1] <- list(seq(10,100,10))

poped.db <- create.poped.database(ff_file="ff",
                                    fg_file="sfg.discrete",
                                    fError_file="feps",
                                    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6), 
                                    sigma=c(0.01,0.25),
                                    groupsize=32,
                                    xt=c( 0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=120,
                                    x=c(70),
                                    discrete_x=x.space)

# use one of the following methods
rs.output <- RS_opt(poped.db,opt_xt=0,opt_x=1,opt_a=0,rsit=rsit)

rs.output <- poped_optimize(poped.db,opt_xt=0,opt_a=0,opt_x=1,rsit=rsit,
                            bUseRandomSearch= 1,
                            bUseStochasticGradient = 0,
                            bUseBFGSMinimizer = 0,
                            bUseLineSearch = 0)

ls.output <- poped_optimize(poped.db,opt_xt=0,opt_a=0, opt_x=1,
                            bUseRandomSearch= 0,
                            bUseStochasticGradient = 0,
                            bUseBFGSMinimizer = 0,
                            bUseLineSearch = 1,
                            ls_step_size=ls_step_size)

mfea.output <- poped_optimize(poped.db,opt_xt=0,opt_a=0,opt_x=1,
                              bUseExchangeAlgorithm=1)

get_rse(mfea.output$fmf,mfea.output$poped.db)


