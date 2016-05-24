## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error to 
##   avoid sample times at very low concentrations (time 0 or very late samoples).
library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.sd.CL

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

######################
# Normal distribution
######################
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed_n <- cbind(ones(length(bpop_vals),1)*1, # normal distribution
                        bpop_vals,
                        ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_n["Favail",]  <- c(0,1,0)
bpop_vals_ed_n

## -- Define initial design  and design space
poped.db.n <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                    fg_file="sfg",
                                    fError_file="feps.add.prop",
                                    bpop=bpop_vals_ed_n, 
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


## ED evaluate using LaPlace approximation 
tic()
output <- evaluate.e.ofv.fim(poped.db.n,use_laplace=TRUE)
toc()
output$E_ofv

\dontrun{
  
  
  ## expected value (roughly)
  tic()
  e_ofv_mc_n <- evaluate.e.ofv.fim(poped.db.n,ED_samp_size=500)
  toc()
  e_ofv_mc_n$E_ofv
  
  
  ## Using ed_laplce_ofv directly
  ed_laplace_ofv(model_switch=poped.db.n$global_model_switch,
                 groupsize=poped.db.n$groupsize,
                 ni=poped.db.n$gni,
                 xtopto=poped.db.n$gxt,
                 xopto=poped.db.n$gx,
                 aopto=poped.db.n$ga,
                 bpopdescr=poped.db.n$gbpop,
                 ddescr=poped.db.n$gd,
                 covd=poped.db.n$covd,
                 sigma=poped.db.n$sigma,
                 docc=poped.db.n$docc, 
                 poped.db.n)
  
  
  ######################
  # Log-normal distribution
  ######################
  
  # Adding 10% log-normal Uncertainty to fixed effects (not Favail)
  bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
  bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                           bpop_vals,
                           ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
  bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
  bpop_vals_ed_ln
  
  ## -- Define initial design  and design space
  poped.db.ln <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                       fg_file="sfg",
                                       fError_file="feps.add.prop",
                                       bpop=bpop_vals_ed_ln, 
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
  
  
  
  ## ED evaluate using LaPlace approximation 
  tic()
  output <- evaluate.e.ofv.fim(poped.db.ln,use_laplace=TRUE)
  toc()
  output$E_ofv
  
  ## expected value (roughly)
  tic()
  e_ofv_mc_ln <- evaluate.e.ofv.fim(poped.db.ln,ED_samp_size=500)[["E_ofv"]]
  toc()
  e_ofv_mc_ln
  
  ## Using ed_laplce_ofv directly
  ed_laplace_ofv(model_switch=poped.db.ln$global_model_switch,
                 groupsize=poped.db.ln$groupsize,
                 ni=poped.db.ln$gni,
                 xtopto=poped.db.ln$gxt,
                 xopto=poped.db.ln$gx,
                 aopto=poped.db.ln$ga,
                 bpopdescr=poped.db.ln$gbpop,
                 ddescr=poped.db.ln$gd,
                 covd=poped.db.ln$covd,
                 sigma=poped.db.ln$sigma,
                 docc=poped.db.ln$docc, 
                 poped.db.ln)
  
  
  
  
}

