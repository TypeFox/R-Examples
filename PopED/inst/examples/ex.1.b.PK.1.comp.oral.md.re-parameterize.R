## using libary models and reparameterizing the problen to KA, KE and V 
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

## -- names match parameters in function defined in ff_file
fg.PK.1.comp.oral.md.param.2 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]), 
                KE=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}

## -- Define design and design space
poped.db <- create.poped.database(ff_file="ff.PK.1.comp.oral.md.KE",
                                    fg_file="fg.PK.1.comp.oral.md.param.2",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    m=2,
                                    sigma=c(0.04,5e-6),
                                    bpop=c(V=72.8,KA=0.25,KE=3.75/72.8,Favail=0.9), 
                                    d=c(V=0.09,KA=0.09,KE=0.25^2), 
                                    notfixed_bpop=c(1,1,1,0),
                                    notfixed_sigma=c(0,0),
                                    xt=c( 1,2,8,240,245),
                                    minxt=c(0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248),
                                    a=cbind(c(20,40),c(24,24)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200,40),
                                    mina=c(0,2))


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
