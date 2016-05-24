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


##-- Model: One comp first order absorption
## -- Analytic solution for both mutiple and single dosing
ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    N = floor(xt/TAU)+1
    y=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
      (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    return(list( y=y,poped.db=poped.db))
  })
}

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}

## -- Residual unexplained variablity (RUV) function
## -- Additive + Proportional  
feps <- function(model_switch,xt,parameters,epsi,poped.db){
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list( y= y,poped.db =poped.db )) 
}

## -- Define design and design space
poped.db <- create.poped.database(ff_file="ff",
                                  fg_file="sfg",
                                  fError_file="feps",
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(0,0),
                                  m=2,
                                  groupsize=20,
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  bUseGrouped_xt=1,
                                  a=cbind(c(20,40),c(24,24)),
                                  maxa=c(200,24),
                                  mina=c(0,24))

##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability 
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of sample times
output <- poped_optimize(poped.db,opt_xt=T,
                         rsit=rsit,sgit=sgit,ls_step_size=ls_step_size,
                         iter_max=iter_max)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db)

# RS+SG+LS optimization of sample times and doses
output <- poped_optimize(poped.db,opt_xt=T,opt_a=T,
                         rsit=rsit,sgit=sgit,ls_step_size=ls_step_size,
                         iter_max=iter_max)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db)

# MFEA optimization with only integers (or multiples of 40 if fast=TRUE) in xt allowed (or original design)
# faster optimization than RS+SG+LS in this case
mfea.output <- poped_optimize(poped.db,opt_xt=T,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=EAStepSize)
get_rse(mfea.output$fmf,mfea.output$poped.db)
plot_model_prediction(mfea.output$poped.db)

# Efficiency of sampling windows
plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1,iNumSimulations=iNumSimulations)
