
library(PopED)

## find the parameters that are needed to define from the structural model
ff.PKPD.1.comp.oral.md.CL.imax
ff.PK.1.comp.oral.md.CL

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU = a[2],
                E0=bpop[5]*exp(b[4]),
                IMAX=bpop[6],
                IC50=bpop[7])
  return( parameters ) 
}



feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  MS <- model_switch
  
  pk.dv <- y*(1+epsi[,1])+epsi[,2]
  pd.dv <-  y*(1+epsi[,3])+epsi[,4]
  
  y[MS==1] = pk.dv[MS==1]
  y[MS==2] = pd.dv[MS==2]
  
  return(list( y= y,poped.db =poped.db )) 
}


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff.PKPD.1.comp.oral.md.CL.imax",
                                  fError_file="feps",
                                  fg_file="sfg",
                                  groupsize=20,
                                  m=3,
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,
                                         E0=1120,IMAX=0.807,IC50=0.0993),  
                                  notfixed_bpop=c(1,1,1,0,1,1,1),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2,E0=0.09), 
                                  sigma=c(0.04,5e-6,0.09,100),
                                  notfixed_sigma=c(0,0,0,0),
                                  xt=c( 1,2,8,240,240,1,2,8,240,240),
                                  minxt=c(0,0,0,240,240,0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248,10,10,10,248,248),
                                  G_xt=c(1,2,3,4,5,1,2,3,4,5),
                                  model_switch=c(1,1,1,1,1,2,2,2,2,2),
                                  a=cbind(c(20,40,0),c(24,24,24)),
                                  bUseGrouped_xt=1,
                                  ourzero=0,
                                  maxa=c(200,40),
                                  mina=c(0,2))


##  create plot of model without variability 
plot_model_prediction(poped.db,facet_scales="free")

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

