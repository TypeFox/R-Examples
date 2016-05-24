library(PopED)
library(deSolve)

# This option is used to make this script run fast but without convergence 
# (fast means a few seconds for each argument at the most).
# This allows you to "source" this file and easily see how things work
# without waiting for more than 10-30 seconds.
# Change to FALSE if you want to run each function so that
# the solutions have converged (can take many minutes).
fast <- TRUE 

# compile and load the qss_one_target.c code.
# to set this up see the 
# "R Package deSolve, Writing Code in Compiled Languages" 
# vingette in the deSolve documentation
system("R CMD SHLIB tmdd_qss_one_target.c")
dyn.load(paste("tmdd_qss_one_target", .Platform$dynlib.ext, sep = ""))

sfg <- function(x,a,bpop,b,bocc){
  parameters=c( CL=bpop[1]*exp(b[1])  ,
                V1=bpop[2]*exp(b[2])  ,
                Q=bpop[3]*exp(b[3])	,
                V2=bpop[4]*exp(b[4])	,
                FAVAIL=bpop[5]*exp(b[5])	,
                KA=bpop[6]*exp(b[6])	,                       
                VMAX=bpop[7]*exp(b[7])	,
                KMSS=bpop[8]*exp(b[8])	,
                R0=bpop[9]*exp(b[9])	,
                KSSS=bpop[10]*exp(b[10])	,
                KDEG=bpop[11]*exp(b[11])	,
                KINT=bpop[12]*exp(b[12])	,
                DOSE=a[1]	,
                SC_FLAG=x[1])   
  return(parameters) 
}

tmdd_qss_one_target_model_compiled <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    
    #The initialization vector for the compartment 
    A_ini <- c(A1=DOSE*SC_FLAG,
               A2=DOSE*(1-SC_FLAG),
               A3=0,
               A4=R0)
    
    #Set up time points for the ODE
    times_xt <- drop(xt)
    times <- sort(times_xt) 
    times <- c(0,times) ## add extra time for start of integration
    
    # solve the ODE    
    out <- ode(A_ini, times, func = "derivs", parms = parameters,
               jacfunc = "jac", # not really needed, speed up is minimal if this is defined or not.
               dllname = "tmdd_qss_one_target",
               initfunc = "initmod", nout = 1, outnames = "Sum")
    
    # extract the time points of the observations
    out = out[match(times_xt,out[,"time"]),]
    
    # Match ODE output to measurements
    RTOT = out[,"A4"]
    CTOT = out[,"A2"]/V1
    CFREE = 0.5*((CTOT-RTOT-KSSS)+sqrt((CTOT-RTOT-KSSS)^2+4*KSSS*CTOT))
    COMPLEX=((RTOT*CFREE)/(KSSS+CFREE))
    RFREE= RTOT-COMPLEX
    
    y[model_switch==1]= RTOT[model_switch==1]
    y[model_switch==2] =CFREE[model_switch==2]
    #y[model_switch==3]=RFREE[model_switch==3]
    
    return(list( y=y,poped.db=poped.db))
  })
}

tmdd_qss_one_target_model_ruv <- function(model_switch,xt,parameters,epsi,poped.db){
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  y[model_switch==1] = log(y[model_switch==1])+epsi[,1]
  y[model_switch==2] = log(y[model_switch==2])+epsi[,2]
  #y[model_switch==3] = log(y[model_switch==3])+epsi[,3]
  
  return(list(y=y,poped.db=poped.db)) 
}

#################################################
# for study 1 in gibiansky,JPKPD,2012 table 2 
#################################################

x.space <- cell(4,1)
for(i in 1:4) x.space[i,1] <- list(c(0,1))


poped.db.1 <- create.poped.database(ff_file="tmdd_qss_one_target_model_compiled",
                                    fError_file="tmdd_qss_one_target_model_ruv",
                                    fg_file="sfg",
                                    groupsize=6,
                                    m=4,      #number of groups
                                    sigma=c(0.04,0.0225), 
                                    bpop=c(CL=0.3,V1=3,Q=0.2,V2=3,FAVAIL=0.7,KA=0.5,VMAX=0,
                                           KMSS=0,R0=0.1,KSSS=0.015,KDEG=10,KINT=0.05),
                                    d=c(CL=0.09,V1=0.09,Q=0.04,V2=0.04,FAVAIL=0.04,KA=0.16,VMAX=0,
                                        KMSS=0,R0=0.09,KSSS=0.09,KDEG=0.04,KINT=0.04), 
                                    notfixed_bpop=c( 1,1,1,1,1,1,0,0,1,1,1,1),
                                    notfixed_d=c( 1,1,1,1,1,1,0,0,1,1,1,1),
                                    xt=c(0.0417,0.25,0.5,1,3,7,14,21,28,35,42,49,56,
                                         0.0417,0.25,0.5,1,3,7,14,21,28,35,42,49,56),
                                    model_switch=c(1,1,1,1,1,1,1,1,1,1,1,1,1,
                                                   2,2,2,2,2,2,2,2,2,2,2,2,2),
                                    bUseGrouped_xt=1,
                                    G_xt=c(1,2,3,4,5,6,7,8,9,10,11,12,13,
                                           1,2,3,4,5,6,7,8,9,10,11,12,13),
                                    a=rbind(100,300,600,1000),
                                    x=rbind(0,0,0,1),
                                    discrete_x=x.space)

plot_model_prediction(poped.db.1,facet_scales="free")

# evaluation time is roughly 3 seconds 
# (macbook pro,OS X 10.10, 2.7 GHz Intel Core i7, 16 GB 1600 MHz DDR3)
# makes optimization reasonable
tic()
FIM <- evaluate.fim(poped.db.1) 
toc() # time is 2.72 sec.
FIM
det(FIM)
get_rse(FIM,poped.db.1)  


#################################################
# for study 1 + 2 in gibiansky,JPKPD,2012 table 2 
#################################################

xt <- zeros(6,30)
study_1_xt <- matrix(rep(c(0.0417,0.25,0.5,1,3,7,14,21,28,35,42,49,56),8),nrow=4,byrow=T)
study_2_xt <- matrix(rep(c(0.0417,1,1,7,14,21,28,56,63,70,77,84,91,98,105),4),nrow=2,byrow=T)
xt[1:4,1:26] <- study_1_xt
xt[5:6,] <- study_2_xt


model_switch <- zeros(6,30)
model_switch[1:4,1:13] <- 1
model_switch[1:4,14:26] <- 2
model_switch[5:6,1:15] <- 1
model_switch[5:6,16:30] <- 2

G_xt <- zeros(6,30)
study_1_G_xt <- matrix(rep(c(1:13),8),nrow=4,byrow=T)
study_2_G_xt <- matrix(rep(c(14:28),4),nrow=2,byrow=T)
G_xt[1:4,1:26] <- study_1_G_xt
G_xt[5:6,] <- study_2_G_xt

x.space.2 <- cell(6,1)
for(i in 1:6) x.space.2[i,1] <- list(c(0,1))

poped.db.2 <- create.poped.database(ff_file="tmdd_qss_one_target_model_compiled",
                                    fError_file="tmdd_qss_one_target_model_ruv",
                                    fg_file="sfg",
                                    groupsize=rbind(6,6,6,6,100,100),
                                    m=6,      #number of groups
                                    sigma=c(0.04,0.0225), 
                                    bpop=c(CL=0.3,V1=3,Q=0.2,V2=3,FAVAIL=0.7,KA=0.5,VMAX=0,
                                           KMSS=0,R0=0.1,KSSS=0.015,KDEG=10,KINT=0.05),
                                    d=c(CL=0.09,V1=0.09,Q=0.04,V2=0.04,FAVAIL=0.04,KA=0.16,VMAX=0,
                                        KMSS=0,R0=0.09,KSSS=0.09,KDEG=0.04,KINT=0.04), 
                                    notfixed_bpop=c( 1,1,1,1,1,1,0,0,1,1,1,1),
                                    notfixed_d=c( 1,1,1,1,1,1,0,0,1,1,1,1),
                                    xt=xt,
                                    model_switch=model_switch,
                                    ni=rbind(26,26,26,26,30,30),
                                    bUseGrouped_xt=1,
                                    G_xt=G_xt,
                                    a=rbind(100,300,600,1000,600,1000),
                                    maxa=1000,
                                    mina=0.1,
                                    x=rbind(0,0,0,1,0,1),
                                    discrete_x=x.space.2)

plot_model_prediction(poped.db.2,facet_scales="free")

# evaluation time is roughly 5 seconds 
# (macbook pro,OS X 10.10, 2.7 GHz Intel Core i7, 16 GB 1600 MHz DDR3)
# makes optimization reasonable
tic()
FIM <- evaluate.fim(poped.db.2) 
toc() 
FIM
det(FIM)
get_rse(FIM,poped.db.2)  # same as in paper: table 1, model 1 
# (except for sigma which appears in the paper to be in SD units, as PFIM reports).

if(!fast){
  # optimize using line search (takes about 2.5 hours)
  ls.output <- poped_optimize(poped.db,opt_xt=0,opt_a=1,
                              bUseRandomSearch= 0,bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,bUseLineSearch = 1,
                              ls_step_size=100)
  
  # Optimized Covariates:
  #   Group 1: 310.069
  # Group 2: 310.069
  # Group 3: 310.069
  # Group 4: 510.049
  # Group 5: 150.085
  # Group 6: 20.098
  # 
  # Efficiency (final_design/initial_design): 1.38869
  # 
  # Expected parameter variance 
  # and relative standard error (%RSE):
  #   Parameter    Values   Variance     RSE   RSE_initial_design
  # bpop[1]    0.3000   7.17e-05    2.82                 2.42
  # bpop[2]    3.0000   3.89e-03    2.08                 2.49
  # bpop[3]    0.2000   2.44e-05    2.47                 2.39
  # bpop[4]    3.0000   2.57e-03    1.69                 2.76
  # bpop[5]    0.7000   2.19e-04    2.11                 3.02
  # bpop[6]    0.5000   4.17e-04    4.09                 4.84
  # bpop[9]    0.1000   4.53e-06    2.13                 2.80
  # bpop[10]    0.0150   1.19e-07    2.30                 3.12
  # bpop[11]   10.0000   1.94e-02    1.39                 2.71
  # bpop[12]    0.0500   7.40e-07    1.72                 2.43
  # D[1,1]    0.0900   1.56e-04   13.89                10.83
  # D[2,2]    0.0900   8.23e-05   10.08                12.27
  # D[3,3]    0.0400   9.49e-05   24.36                22.27
  # D[4,4]    0.0400   2.85e-05   13.36                20.03
  # D[5,5]    0.0400   4.11e-05   16.03                24.24
  # D[6,6]    0.1600   5.76e-04   15.00                19.41
  # D[9,9]    0.0900   8.79e-05   10.42                11.89
  # D[10,10]    0.0900   1.18e-04   12.09                13.45
  # D[11,11]    0.0400   1.63e-05   10.09                20.40
  # D[12,12]    0.0400   2.76e-05   13.14                18.13
  # SIGMA[1,1]    0.0400   1.17e-06    2.70                 2.73
  # SIGMA[2,2]    0.0225   4.89e-07    3.11                 2.91
  
  # Total running time: 10767.9 seconds
  
  plot_model_prediction(ls.output$poped.db,facet_scales="free")
}
