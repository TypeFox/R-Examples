library(PopED)

# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 

rsit <- ifelse(fast,3,300)
ED_samp_size <- ifelse(fast,20,100)


ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption + inhibitory imax
  ## -- works for both mutiple and single dosing  
  with(as.list(parameters),{
    
    y=xt
    MS <- model_switch
    
    # PK model
    N = floor(xt/TAU)+1
    CONC=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
      (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    
    # PD model
    EFF = E0*(1 - CONC*IMAX/(IC50 + CONC))
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

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
  returnArgs <- ff(model_switch,xt,parameters,poped.db) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  MS <- model_switch
  
  pk.dv <- y*(1+epsi[,1])+epsi[,2]
  pd.dv <-  y*(1+epsi[,3])+epsi[,4]
  
  y[MS==1] = pk.dv[MS==1]
  y[MS==2] = pd.dv[MS==2]
  
  return(list( y= y,poped.db =poped.db )) 
}

# Adding 10% Uncertainty to IC50 parameter
bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,E0=1120,IMAX=0.807,IC50=0.0993)
bpop_vals_ed <- cbind(zeros(7,1),bpop_vals,zeros(7,1)) 
bpop_vals_ed["IC50",1] <- 1 # normal distrtibution
bpop_vals_ed["IC50",3] <- (bpop_vals_ed["IC50",2]*0.1)^2
bpop_vals_ed

poped.db <- create.poped.database(ff_file="ff",
                                  fError_file="feps",
                                  fg_file="sfg",
                                  groupsize=20,
                                  m=3,
                                  bpop=bpop_vals_ed,  
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



## ED evaluate using API. 
output <- evaluate.e.ofv.fim(poped.db,ofv_calc_type=4,ED_samp_size = ED_samp_size)
output$E_ofv


## optimization with random search
output <- poped_optimize(poped.db,opt_xt=T, opt_a=T,
                         d_switch=F,ED_samp_size=ED_samp_size,
                         rsit=rsit,
                         ofv_calc_type=4)


get_rse(output$fmf,output$poped.db)
result.db <- output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T,separate.groups=T,facet_scales="free")






