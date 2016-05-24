library(PopED)
library(deSolve)

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


PK.1.comp.oral.sd.fg.param.1 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}

bpop.vals.param.1 <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
d.vals.param.1 <- c(V=0.09,KA=0.09,CL=0.25^2)
notfixed_bpop <- c(1,1,1,0)

PK.1.comp.oral.sd.fg.param.2 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]), 
                KE=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  return( parameters ) 
}

bpop.vals.param.2 <- c(V=72.8,KA=0.25,KE=3.75/72.8,Favail=0.9)
d.vals.param.2 <- c(V=0.09,KA=0.09,KE=0.25^2)
notfixed_bpop <- c(1,1,1,0)


PK.1.comp.oral.sd.discrete.fg.param.1 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=x[1])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}


PK.1.comp.oral.sd.ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
    return(list( y= y,poped.db=poped.db))
  })
}

PK.1.comp.oral.sd.ff.ode <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE, A2 = 0)
    times <- drop(xt)
    times <- sort(times) 
    times <- c(0,times) # add extra time for start of integration
    out   <- ode(A_ini, times, PK.1.comp.oral.ode, parameters) #,atol=1e-13,rtol=1e-13)
    y = out[,"A2"]/(V/Favail)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) ## must be a row vector
    return(list( y= y,poped.db=poped.db)) 
  })
}

PK.1.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dA1  <- -KA*A1
    dA2  <- KA*A1 - KE*A2
    
    return(list(c(dA1, dA2)))
  })
}

poped.db.1 <- create.poped.database(ff_file="PK.1.comp.oral.sd.ff",
                                    fg_file="PK.1.comp.oral.sd.fg.param.1",
                                    fError_file="feps.add.prop",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,CL=0.15,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=notfixed_bpop,
                                    xt=c(0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))

poped.db.2 <- create.poped.database(ff_file="PK.1.comp.oral.sd.ff",
                                    fg_file="PK.1.comp.oral.sd.fg.param.2",
                                    fError_file="feps.add.prop",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,KE=0.15/8,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=notfixed_bpop,
                                    xt=c(0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))

poped.db.3 <- create.poped.database(ff_file="PK.1.comp.oral.sd.ff.ode",
                                    fError_file="feps.add.prop",
                                    fg_file="PK.1.comp.oral.sd.fg.param.1",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,CL=0.15,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=notfixed_bpop,
                                    xt=c(0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))

for(i in 1:3){
  poped.db <- eval(parse(text=paste("poped.db.",i,sep="")))
  
  ##  create plot of model 
  print(plot_model_prediction(poped.db,IPRED=T,DV=T))
  
  ## evaluate initial design
  FIM <- evaluate.fim(poped.db) 
  print(FIM)
  print(det(FIM))
  print(get_rse(FIM,poped.db))
  
  if(i<3){
    # RS+SG+LS optimization of sample times
    output <- poped_optimize(poped.db,opt_xt=T,
                             rsit=rsit,sgit=sgit,ls_step_size=ls_step_size,
                             iter_max=iter_max)
    print(get_rse(output$fmf,output$poped.db))
    print(plot_model_prediction(output$poped.db))
    
    
    # MFEA optimization with only integers (or multiples of 40 if fast=TRUE) in xt allowed (or original design)
    # faster optimization than RS+SG+LS in this case
    mfea.output <- poped_optimize(poped.db,opt_xt=T,
                                  bUseExchangeAlgorithm=1,
                                  EAStepSize=EAStepSize)
    get_rse(mfea.output$fmf,mfea.output$poped.db)
    plot_model_prediction(mfea.output$poped.db)
    
    # Efficiency of sampling windows
    plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1,iNumSimulations=iNumSimulations)
  }
}
