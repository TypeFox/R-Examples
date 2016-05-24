## HCV example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)
library(deSolve)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c(p=bpop[1],
               d=bpop[2],
               e=bpop[3],
               s=bpop[4],
               KA=bpop[5] + b[1],
               KE=bpop[6] + b[2],
               VD=bpop[7] + b[3],
               EC50=bpop[8] + b[4],
               n=bpop[9] + b[5],
               delta=bpop[10] + b[6],
               c=bpop[11] + b[7],
               DOSE=a[1],
               TINF=a[2],
               TAU=a[3])
  return(parameters)
}

#' To make evaluation time more reasonable we use compiled code
#' To set this up see the 
#' "R Package deSolve, Writing Code in Compiled Languages" 
#' vingette in the deSolve documentation
#' make sure your working directory is where this file is located
system("R CMD SHLIB HCV_ode.c")
dyn.load(paste("HCV_ode", .Platform$dynlib.ext, sep = ""))

ff_ODE_compiled <- function(model_switch,xt,parameters,poped.db){
  parameters[5:11] <- exp(parameters[5:11])
  with(as.list(parameters),{
    A_ini  <- c(A1 = 0, A2 = 0, A3=c*delta/(p*e), 
                A4=(s*e*p-d*c*delta)/(p*delta*e),
                A5=(s*e*p-d*c*delta)/(c*delta*e))
    
    #Set up time points for the ODE
    times_xt <- drop(xt)
    times <- c(0,times_xt) ## add extra time for start of integration
    times <- sort(times) 
    times <- unique(times) # remove duplicates
    
    # compute values from ODEs
    out <- ode(A_ini, times, func = "derivs", parms = parameters, 
               #method="daspk",
               #jacfunc = "jac", 
               dllname = "HCV_ode",
               initfunc = "initmod", #nout = 1, outnames = "Sum",
               atol=1e-12,rtol=1e-12)
    
    # grab timepoint values
    out = out[match(times_xt,out[,"time"]),]
    
    y <- xt*0
    pk <- out[,"A2"]/VD
    pd <- out[,"A5"]
    y[model_switch==1] <- pk[model_switch==1]
    y[model_switch==2] <- log10(pd[model_switch==2])
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}

feps_ODE_compiled <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  MS<-model_switch
  y <- ff_ODE_compiled(model_switch,xt,parameters,poped.db)[[1]]
  
  pk.dv <- y + epsi[,1]
  pd.dv <- y + epsi[,2]
  
  y[MS==1] = pk.dv[MS==1]
  y[MS==2] = pd.dv[MS==2]
  
  return(list(y=y,poped.db=poped.db))
}

## -- Define initial design  and design space
poped_db_compiled <- create.poped.database(ff_file="ff_ODE_compiled",
                                           fg_file="sfg",
                                           fError_file="feps_ODE_compiled",
                                           bpop=c(p=100,
                                                  d=0.001,
                                                  e=1e-7,
                                                  s=20000,
                                                  KA=log(0.8),
                                                  KE=log(0.15),
                                                  VD=log(100), 
                                                  #VD=log(100000),
                                                  EC50=log(0.12), 
                                                  #EC50=log(0.00012),
                                                  n=log(2),
                                                  delta=log(0.2),
                                                  c=log(7)),
                                           notfixed_bpop=c(0,0,0,0,1,1,1,1,1,1,1),
                                           d=c(KA=0.25,
                                               KE=0.25,
                                               VD=0.25,
                                               EC50=0.25,
                                               n=0.25,
                                               delta=0.25,
                                               c=0.25),
                                           sigma=c(0.04,0.04),
                                           groupsize=30,
                                           xt=c(0,0.25,0.5,1,2,3,4,7,10,14,21,28,
                                                0,0.25,0.5,1,2,3,4,7,10,14,21,28),
                                           model_switch=c(rep(1,12),rep(2,12)),
                                           a=c(180,1,7))

##  create plot of model without variability 
plot_model_prediction(poped_db_compiled, facet_scales = "free")

#' #####################################
#' The reduced FIM 
#' ####################################
tic() # computation time
FIM_compiled <- evaluate.fim(poped_db_compiled) 
toc()

#' design evaluation
crit <- det(FIM_compiled)^(1/length(get_unfixed_params(poped_db_compiled)[["all"]]))
crit
rse <- get_rse(FIM_compiled,poped_db_compiled) # this is for the log of the fixed effect parameters
rse_norm <- sqrt(diag(inv(FIM_compiled)))*100 # this is approximately the RSE for the normal scale of the fixed effects 
rse[1:7] <- rse_norm[1:7] # replace the log scale for the normal scale RSE values
rse

#' Evaluation comparison to 
#' Nyberg et al., "Methods and software tools for design evaluation 
#' for population pharmacokinetics-pharmacodynamics studies", 
#' Br. J. Clin. Pharm., 2014. 
crit_reference_reduced <- 248.8
rse_reference_reduced <- c(12.1,10.5,10.0,15.8,10.4,9.4,11.0,40.0,30.8,28.8,60.4,28.8,27.2,32.8,8.5,9.0)

# the relative differences in percent
(rse - rse_reference_reduced)/rse_reference_reduced * 100
(crit - crit_reference_reduced)/crit_reference_reduced * 100

#' #####################################
#' The full FIM. 
#' #####################################
FIM_compiled_full <- evaluate.fim(poped_db_compiled,fim.calc.type = 0) 

#' design evaluation
crit_full <- det(FIM_compiled_full)^(1/length(get_unfixed_params(poped_db_compiled)[["all"]]))
crit_full
rse_full <- get_rse(FIM_compiled_full,poped_db_compiled) # this is for the log of the fixed effect parameters
rse_norm_full <- sqrt(diag(inv(FIM_compiled_full)))*100 # this is approximately the RSE for the normal scale of the fixed effects 
rse_full[1:7] <- rse_norm_full[1:7] # replace the log scale for the normal scale RSE values
rse_full

#' Evaluation compared to 
#' Nyberg et al., "Methods and software tools for design evaluation 
#' for population pharmacokinetics-pharmacodynamics studies", 
#' Br. J. Clin. Pharm., 2014. 
crit_reference_full <- 318.2
rse_reference_full <- c(8.6,6.9,8.4,13.5,7.5,8.5,8.7,43.2,37.2,33.2,66.4,32.8,31.6,33.6,8.5,9.3)

# the relative differences in percent
(rse_full - rse_reference_full)/rse_reference_full * 100
(crit_full - crit_reference_full)/crit_reference_full * 100



