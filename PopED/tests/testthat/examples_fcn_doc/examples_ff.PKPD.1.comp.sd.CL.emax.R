
library(PopED)

## find the parameters that are needed to define from the structural model
ff.PKPD.1.comp.sd.CL.emax

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( 
    CL=bpop[1]*exp(b[1])  ,
    V=bpop[2]*exp(b[2])  ,
    E0=bpop[3]*exp(b[3])  ,
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

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff.PKPD.1.comp.sd.CL.emax",
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


##  create plot of model without variability 
plot_model_prediction(poped.db,facet_scales="free")

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

