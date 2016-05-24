library(PopED)

library(deSolve)

# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 


PK.1.comp.oral.md.ff.ode <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt) #xt[,,drop=T] 
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, PK.1.comp.oral.ode, parameters, events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V/Favail)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

PK.1.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {    
    dA1 <- -KA*A1
    dA2 <- KA*A1 - (CL/V)*A2
    return(list(c(dA1, dA2)))
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

poped.db <- create.poped.database(ff_file="PK.1.comp.oral.md.ff.ode",
                                    fError_file="feps",
                                    fg_file="sfg",
                                    groupsize=20,
                                    m=2,      #number of groups
                                  sigma=c(0.04,5e-6),
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  notfixed_bpop=c(1,1,1,0),
                                  notfixed_sigma=c(0,0),
                                    xt=c( 1,2,8,240,245),
                                    minxt=c(0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248),
                                    a=cbind(c(20,40),c(24,24)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200,40),
                                    mina=c(0,2))


plot_model_prediction(poped.db)

FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

if(!fast){
  # optimize using line search
  ls.output <- poped_optimize(poped.db,opt_xt=1,
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 1)
  
  plot_model_prediction(ls.output$poped.db)
  
}