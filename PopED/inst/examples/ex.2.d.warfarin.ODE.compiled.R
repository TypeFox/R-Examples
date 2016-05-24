## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error to 
##   avoid sample times at very low concentrations (time 0 or very late samoples).

## Model described with an ODE
library(PopED)
library(deSolve)

# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 

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

ff.ODE <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE, A2 = 0)
    times <- drop(xt)##xt[,,drop=T] 
    times <- sort(times) 
    times <- c(0,times) ## add extra time for start of integration
    out   <- ode(A_ini, times, one.comp.ode, parameters)#,atol=1e-13,rtol=1e-13)
    y = out[,"A2"]/(V/Favail)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}

one.comp.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {    
    dA1  <- -KA*A1
    dA2  <- KA*A1 - CL/V*A2
    return(list(c(dA1, dA2)))
  })
}


feps.ODE <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional + additive
  y <- ff.ODE(model_switch,xt,parameters,poped.db)[[1]] 
  y = y*(1+epsi[,1]) + epsi[,2]  
  return(list(y=y,poped.db=poped.db)) 
}


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff.ODE",
                                  fg_file="sfg",
                                  fError_file="feps.ODE",
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

# MFEA optimization with only integer times and doses allowed
# speed is quite slow with the ODE model if not using compiled code (see below)
if(!fast){
  mfea.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(mfea.output$fmf,mfea.output$poped.db)
  plot_model_prediction(mfea.output$poped.db)
}


##################################
# Compiled ODE
##################################

# compile and load the qss_one_target.c code.
# to set this up see the 
# "R Package deSolve, Writing Code in Compiled Languages" 
# vingette in the deSolve documentation

system("R CMD SHLIB one_comp_oral_CL.c")
dyn.load(paste("one_comp_oral_CL", .Platform$dynlib.ext, sep = ""))

ff.ODE.compiled <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE, A2 = 0)
    times <- drop(xt)##xt[,,drop=T] 
    times <- sort(times) 
    times <- c(0,times) ## add extra time for start of integration
    out <- ode(A_ini, times, func = "derivs", parms = parameters,
               #jacfunc = "jac", # not really needed, speed up is minimal if this is defined or not.
               dllname = "one_comp_oral_CL",
               initfunc = "initmod", nout = 1, outnames = "Sum")    
    y = out[,"A2"]/(V/Favail)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}

feps.ODE.compiled <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional + additive
  y <- ff.ODE.compiled(model_switch,xt,parameters,poped.db)[[1]] 
  y = y*(1+epsi[,1]) + epsi[,2]  
  return(list(y=y,poped.db=poped.db)) 
}

## -- Define initial design  and design space
poped.db.compiled <- create.poped.database(ff_file="ff.ODE.compiled",
                                           fg_file="sfg",
                                           fError_file="feps.ODE.compiled",
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
plot_model_prediction(poped.db.compiled)

##  create plot of model with variability 
plot_model_prediction(poped.db.compiled,IPRED=T,DV=T)

## evaluate initial design
FIM.compiled <- evaluate.fim(poped.db.compiled) 
det(FIM.compiled)
get_rse(FIM.compiled,poped.db.compiled)

# no difference in computation
det(FIM)-det(FIM.compiled)

# but a huge different in computation time (22 times faster with the compiled code)
if(!fast){
  library(microbenchmark)
  compare <- microbenchmark(evaluate.fim(poped.db.compiled), evaluate.fim(poped.db), times = 100)
  compare
  library(ggplot2)
  autoplot(compare)
}

## making optimization times resonable
output <- poped_optimize(poped.db.compiled,opt_xt=T, opt_a=T,
                         rsit=rsit,sgit=sgit,ls_step_size=ls_step_size,
                         iter_max=iter_max)


##################################
# comapre to the analytic solution
##################################
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


poped.db.1 <- create.poped.database(ff_file="ff",
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

FIM.1 <- evaluate.fim(poped.db.1) 
FIM.1
det(FIM.1)
get_rse(FIM.1,poped.db.1)

## differences can be reduced by decreasing atol and rtol in the "ode" function 
(det(FIM) - det(FIM.1))/det(FIM)*100
(det(FIM)^(1/size(FIM,1)) - det(FIM.1)^(1/size(FIM.1,1)))/det(FIM)^(1/size(FIM,1))*100
(FIM - FIM.1)/FIM*100


## computation times are 6x faster with the analytic solution
if(!fast){
  library(microbenchmark)
  compare <- microbenchmark(evaluate.fim(poped.db.compiled), evaluate.fim(poped.db.1), times = 100)
  compare
  library(ggplot2)
  autoplot(compare)
}
