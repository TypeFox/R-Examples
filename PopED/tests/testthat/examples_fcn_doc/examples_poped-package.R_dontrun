
library(PopED)

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
                                  groupsize=20,
                                  m=2,
                                  sigma=diag(c(0.04,5e-6)),
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  notfixed_bpop=c(1,1,1,0),
                                  notfixed_sigma=c(0,0),
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  a=cbind(c(20,40),c(24,24)),
                                  bUseGrouped_xt=1,
                                  maxa=c(200,24),
                                  mina=c(0,24))

##  create plot of model without variability 
plot_model_prediction(poped.db)


\dontrun{
  
  ##  create plot of model with variability 
  plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)

}

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

\dontrun{
  
  # Optimization of sample times
  output <- poped_optim(poped.db, opt_xt=TRUE, parallel = TRUE)
  get_rse(output$FIM, output$poped.db)
  plot_model_prediction(output$poped.db, IPRED=F, DV=F)
  
  # Optimization of sample times and doses
  output <- poped_optim(poped.db, opt_xt=TRUE, opt_a=TRUE, parallel = TRUE)
  get_rse(output$FIM,output$poped.db)
  plot_model_prediction(output$poped.db,IPRED=F,DV=F)
  
  # Optimization with only integer times and doses allowed
  mfea.output <- poped_optimize(poped.db,opt_xt=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(mfea.output$fmf,mfea.output$poped.db)
  plot_model_prediction(mfea.output$poped.db)
  
  # Efficiency of sampling windows
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5)
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1)
  
}
