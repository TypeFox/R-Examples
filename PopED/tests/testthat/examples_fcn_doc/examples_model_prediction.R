## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.md.CL

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

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                  fg_file="sfg",
                                  fError_file="feps.prop",
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=0.01,
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70)

## data frame with model predictions
df_1 <- model_prediction(poped.db)

##  data frame with variability 
df_2 <- model_prediction(poped.db,DV=TRUE)

## data frame with variability (only IPRED, no DV)
df_3 <- model_prediction(poped.db,IPRED=TRUE)

## data frame with model predictions, no continuous design variables in data frame
df_4 <- model_prediction(poped.db,include_a = FALSE)

## -- 2 groups
poped.db.2 <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                    fg_file="sfg",
                                    fError_file="feps.prop",
                                    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6), 
                                    sigma=0.01,
                                    groupsize=rbind(3,3),
                                    m=2,
                                    xt=c( 0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=120,
                                    a=rbind(70,50))

df_5 <- model_prediction(poped.db.2,DV=TRUE)

## without a poped database, just describing the design
## Useful for creating datasets for use in other software (like NONMEM)
design_1 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  m=2,
  groupsize=3)

design_2 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  m=2,
  groupsize=3,
  a=c(WT=70,AGE=50))

design_3 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  m=2,
  groupsize=3,
  a=list(c(WT=70,AGE=50),c(AGE=45,WT=60)))

df_6 <- model_prediction(design=design_1)
df_7 <- model_prediction(design=design_2)
df_8 <- model_prediction(design=design_3)
df_9 <- model_prediction(design=design_3,DV=TRUE)

# generate random deviations in WT for each individual
df_10 <- model_prediction(design=design_3,DV=TRUE,
                          manipulation=expression({for(id in unique(ID)) 
                              WT[ID==id] = rnorm(1,WT[ID==id],WT[ID==id]*0.1);id <- NULL}))

# generate random deviations in WT and AGE for each individual
df_11 <- model_prediction(design=design_3,DV=TRUE,
                          manipulation=list(
                            expression(for(id in unique(ID)) 
                              WT[ID==id] = rnorm(1,WT[ID==id],WT[ID==id]*0.1)),
                            expression(for(id in unique(ID)) 
                              AGE[ID==id] = rnorm(1,AGE[ID==id],AGE[ID==id]*0.2)),
                            expression(id <- NULL)
                          ))

## create dosing rows 
dosing_1 <- list(list(AMT=1000,RATE=NA,Time=0.5),list(AMT=3000,RATE=NA,Time=0.5))
dosing_2 <- list(list(AMT=1000,RATE=NA,Time=0.5))
dosing_3 <- list(list(AMT=1000,Time=0.5))
dosing_4 <- list(list(AMT=c(1000,20),Time=c(0.5,10))) # multiple dosing


df_12 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_1)
df_13 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_2)
df_14 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_3)
df_15 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_4)

df_16 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_4,filename="test.csv")

model_prediction(design=design_3,DV=TRUE,dosing=dosing_4,model_num_points = 10)
model_prediction(design=design_3,DV=TRUE,dosing=dosing_4,model_num_points = 10,model_minxt=20)

design_4 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  model_switch=c(1,1,1,1,2,2,2,2),
  m=2,
  groupsize=3,
  a=list(c(WT=70,AGE=50),c(AGE=45,WT=60)))

model_prediction(design=design_4,DV=TRUE,dosing=dosing_4)
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = 10)
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = 10,
                 model_minxt=10,model_maxxt=100)
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = 10,
                 model_minxt=c(20,20),model_maxxt=c(100,100))
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = c(10,10),
                 model_minxt=c(20,20),model_maxxt=c(100,100))
