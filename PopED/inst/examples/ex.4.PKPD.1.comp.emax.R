library(PopED)

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

ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    MS <- model_switch
    
    # PK model
    CONC = DOSE/V*exp(-CL/V*xt) 
    
    # PD model
    EFF = E0 + CONC*EMAX/(EC50 + CONC)
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( 
    CL=bpop[1]*exp(b[1])  ,
    V=bpop[2]*exp(b[2])	,
    E0=bpop[3]*exp(b[3])	,
    EMAX=bpop[4]*exp(b[4])	,
    EC50=bpop[5]*exp(b[5])	,
    DOSE=a[1]
  )
  return( parameters ) 
}

feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional PK + additive PD
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  MS <- model_switch
  
  prop.err <- y*(1+epsi[,1])
  add.err <- y+epsi[,2]
  
  y[MS==1] = prop.err[MS==1]
  y[MS==2] = add.err[MS==2]
  
  return(list( y= y,poped.db =poped.db )) 
}

poped.db <- create.poped.database(ff_file="ff",
                                  fError_file="feps",
                                  fg_file="sfg",
                                  groupsize=20,
                                  m=3,
                                  sigma=diag(c(0.15,0.15)),
                                  bpop=c(CL=0.5,V=0.2,E0=1,EMAX=1,EC50=1),  
                                  d=c(CL=0.01,V=0.01,E0=0.01,EMAX=0.01,EC50=0.01), 
                                  xt=c( 0.33,0.66,0.9,5,0.1,1,2,5),
                                  model_switch=c( 1,1,1,1,2,2,2,2),
                                  minxt=0,
                                  maxxt=5,
                                  a=rbind(2.75,5,10),
                                  bUseGrouped_xt=1,
                                  maxa=10,
                                  mina=0.1)


plot_model_prediction(poped.db,facet_scales="free")
plot_model_prediction(poped.db,IPRED=T,DV=T,facet_scales="free",separate.groups=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of sample times and doses
output <- poped_optimize(poped.db,opt_xt=T,opt_a=T,
                         rsit=rsit,sgit=sgit,ls_step_size=ls_step_size,
                         iter_max=iter_max)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db,facet_scales="free")


