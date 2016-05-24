## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Evaluating with uncertainty around parameter values in the model

library(PopED)

# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 

EAStepSize <- ifelse(fast,60,1)
rsit <- ifelse(fast,3,300)

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

# Adding 10% Uncertainty to all fixed effects (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                      bpop_vals,
                      ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed["Favail",] <- c(0,1,0)
bpop_vals_ed


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff",
                                  fg_file="sfg",
                                  fError_file="feps",
                                  bpop=bpop_vals_ed, 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70,
                                  mina=0,
                                  maxa=100,
                                  ED_samp_size=20)

## ED evaluate.
## result is inaccurate (run several times to see)
## increase ED_samp_size for a more accurate calculation
output <- evaluate.e.ofv.fim(poped.db,ED_samp_size=20)
output$E_ofv
output$E_fim

## optimization with random search 
output <- poped_optimize(poped.db,opt_xt=T, opt_a=T,
                         d_switch=F,ED_samp_size=20,
                         rsit=rsit)

# API optimization: E(ln(det(FIM)))
output <- poped_optimize(poped.db,opt_xt=T, opt_a=T,
                         d_switch=F,ED_samp_size=20,
                         rsit=rsit,
                         ofv_calc_type=4)

## MFEA optimization
mfea.output <- poped_optimize(poped.db,
                              opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=EAStepSize,
                              ED_samp_size=20,
                              d_switch=0)

## ED evaluation of ofv using Laplace approximation 
## deterministic calculation, relatively fast
## can be more stable for optimization
tic()
output <- evaluate.e.ofv.fim(poped.db,use_laplace=TRUE)
toc()
output$E_ofv

## optimization with random search and Laplace
output <- poped_optimize(poped.db,opt_xt=T, opt_a=T,
                         d_switch=F,use_laplace=T,
                         rsit=rsit)

