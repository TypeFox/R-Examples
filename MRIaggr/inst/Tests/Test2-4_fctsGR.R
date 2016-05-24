#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Fichier test 2.4 : Test des fonctions GR 
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require(MRIaggr)

options(error=function() traceback(2)) 
options(max.print=10000)

#### 1- calc ####

#### >a calcCriteriaGR - private #### 
# ?calcCriteriaGR
# MRIaggr:::calcCriteriaGR
# function(contrast,groups,W=NULL,sigma=NULL,breaks,
# criterion_d1=FALSE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE){
  
#### test
data(MRIaggr.Pat1_red, package="MRIaggr")

carto <- selectContrast(MRIaggr.Pat1_red,num=3,param=c("TTP_t0","MASK_DWI_t0"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,num=3,hemisphere="lesion")
W <- calcW(coords,range=sqrt(2))

indexN <- which(carto$MASK_DWI_t0==1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

resGR_0.8 <- calcGR(contrast=carto$TTP_t0, W=W, 
                 seed=seed, 
                 sigma_max=0.8, 
                 iter_max=20)

resEvalGR_0.8 <- MRIaggr:::calcCriteriaGR(contrast=carto$TTP_t0,groups=1:nrow(carto) %in% resGR_0.8$GR,W=W,sigma=0.8,breaks=resGR_0.8$breaks,
                            criterion_d1=TRUE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE)
resEvalGR_0.8

resGR_1 <- calcGR(contrast=carto$TTP_t0, W=W, 
                seed=seed, 
                sigma_max=1, 
                iter_max=20)

resEvalGR_1 <- MRIaggr:::calcCriteriaGR(contrast=carto$TTP_t0,groups=1:nrow(carto) %in% resGR_1$GR,W=W,sigma=1,breaks=resGR_1$breaks,
                                        criterion_d1=TRUE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE)
resEvalGR_1


resGR_2 <- calcGR(contrast=carto$TTP_t0, W=W, 
                  seed=seed, 
                  sigma_max=2, 
                  iter_max=20)

resEvalGR_2 <- MRIaggr:::calcCriteriaGR(contrast=carto$TTP_t0,groups=1:nrow(carto) %in% resGR_2$GR,W=W,sigma=2,breaks=resGR_2$breaks,
                                        criterion_d1=TRUE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE)
resEvalGR_2

#### example

#### >b calcGR #### 
# ?calcGR
# calcGR
# function(contrast,W,seed,sigma_max,range=c(-Inf,+Inf),range.seed=c(-Inf,+Inf),coords=NULL,sd.iter=4,                   
# nbreaks=100,breaks=NULL,scale=FALSE,iter_max=100,sd.robust=FALSE,trace=TRUE,
# history_sigma=FALSE,history_step=FALSE,history_front=FALSE)

#### test
resGR0 <- calcGR(contrast=carto$TTP_t0, W=W, 
                 seed=seed, 
                 sigma_max=1, 
                 iter_max=20,
                 range=c(7,10))

# display
multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,hemisphere="lesion",legend=FALSE,
             breaks=seq(0,10,0.1),as.logical=TRUE,ylim=c(25,45), xlim=c(75,95),cex=2,            
             index1=list(coords=coords[resGR0$GR,],pch=20,cex=1),
             index2=list(coords=coords[seed,],pch=20,cex=1)
)

#### test sphere
# n <- 100
# for(iter in c(10,20,30,40,50)){
#   radius <- iter
#   coords <- cbind(expand.grid(1:n,1:n),1)
#   contrast <- as.numeric(sqrt((coords[,1]-n/2)^2+(coords[2]-n/2)^2)<radius)
#   # W <- calcW(coords,range=sqrt(2),upper=NULL)
#   
#   # multiplot(coords,
#   #              contrast)
#   
#   resGR0 <- calcGR(contrast=contrast, W=W, 
#                    seed=which( (coords[,1]==50)*(coords[,2]==50) == 1 ),
#                    sigma_max=10,range=c(0.5,+Inf),
#                    iter_max=50,history_sigma=FALSE,history_step=TRUE,history_front=TRUE)
#   
#   # multiplot(coords,contrast,
#   #              index1=coords[resGR0$GR,],main=iter)
#   
#   mean_step <- mean(resGR0$history_step)
#   sd_step <- sd(resGR0$history_step)
#   
#   multiplot(coords[resGR0$GR,],(resGR0$history_step)/sd_step,
#                main=iter)
# }
# hist(resGR0$history_step)
# 
# IC_step <- mean_step + c(-1,1)*sd_step*qnorm(0.975)


#### example
## load an MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

## display raw parameter
multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,
             breaks=seq(0,10,0.1),as.logical=TRUE,
             index1=list(coords="MASK_DWI_t0",outline=TRUE))

## extract raw parameter, coordinates and compute the neighborhood matrix
carto <- selectContrast(MRIaggr.Pat1_red,num=3,param=c("TTP_t0","MASK_DWI_t0"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,num=3,hemisphere="lesion")
W <- calcW(coords,range=sqrt(2))

## the seed is taken to be the point with the largest TTP in the lesion mask
indexN <- which(carto$MASK_DWI_t0==1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

## Display step by step the GR algorithm with sigma = 1
for(iter in c(0,1,2,5,10,15)){
  resGR1 <- calcGR(contrast=carto$TTP_t0, W=W, 
                   seed=seed, 
                   sigma_max=1, 
                   iter_max=iter,trace=FALSE)
  
  multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,hemisphere="lesion",legend=FALSE,
               breaks=seq(0,10,0.1),as.logical=TRUE,ylim=c(25,45), xlim=c(75,95),cex=2,            
               index1=list(coords=coords[resGR1$GR,],pch=20,cex=1),
               index2=list(coords=coords[seed,],pch=20,cex=1)
  )
}

## GR with sigma = 0.8
resGR2 <- calcGR(contrast=carto$TTP_t0, W=W, 
                 seed=seed,sigma_max=0.8,iter_max=50,
                 history_step=TRUE,history_front=TRUE)

## display 
# display the GR over the raw contrast
multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,hemisphere="lesion",legend=FALSE,
             breaks=seq(0,10,0.1),as.logical=TRUE,ylim=c(25,45), xlim=c(75,95),cex=2,            
             index1=list(coords=coords[resGR2$GR,],pch=20,cex=1)
)

# display the step of inclusion in GR group for each observation
multiplot(coords[resGR2$GR,],
             resGR2$history_step,breaks=0:10,
             index1=list(coords=coords[seed,]),
             palette=rainbow(10)
)

# display the front propagation 
multiplot(coords[resGR2$GR,],
             resGR2$Mfront[,10],
             index1=list(coords=coords[seed,])
)

## GR with  sigma = 2
resGR3 <- calcGR(contrast=carto$TTP_t0, W=W, 
                 seed=seed, 
                 sigma_max=2, 
                 iter_max=20)

## display
multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,hemisphere="lesion",legend=FALSE,
             breaks=seq(0,10,0.1),as.logical=TRUE,ylim=c(25,45), xlim=c(75,95),cex=2,            
             index1=list(coords=coords[resGR3$GR,],pch=20,cex=1)
)




#### >c calcSigmaGR  #### 
# ?calcSigmaGR
# calcSigmaGR
# function(contrast,W,seed,sigma,range=c(-Inf,+Inf),range.seed=c(-Inf,+Inf),
# breaks=100,scale=FALSE,iter_max=100,sd.robust=FALSE,
# criterion_d1=FALSE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE,
# trace=TRUE,mar=rep(2,4),mgp=c(2,0.75,0),main="",
# window=FALSE,filename="calcSigmaGR",width=1000,height=700,path=NULL,unit="px",res=NA)

#### test

#### example
## load an MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

## select data
carto <- selectContrast(MRIaggr.Pat1_red,num=3,param=c("TTP_t0","MASK_DWI_t0"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,num=3,hemisphere="lesion")
W <- calcW(coords,range=sqrt(2))

indexN <- which(carto$MASK_DWI_t0==1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

## find optimal sigma
resGR_sigma <- calcSigmaGR(contrast=carto$TTP_t0, W=W,seed=seed,
                           sigma=seq(0.8,2,0.1),iter_max=50)

  
#### >d GRalgo - private #### 
# ?GRalgo
# MRIaggr:::GRalgo
# function (contrast, W, seed, sigma_max, range, breaks, 
# step, operator, iter_max, history)  

#### test
data(MRIaggr.Pat1_red, package="MRIaggr")

multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,
             breaks=seq(0,10,0.1),as.logical=TRUE,
             index1=list(coords="MASK_DWI_t0",outline=TRUE))

carto <- selectContrast(MRIaggr.Pat1_red,num=3,param=c("TTP_t0","MASK_DWI_t0","T2_FLAIR_t2"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,num=3,hemisphere="lesion")
W <- calcW(coords,range=sqrt(2))

indexN <- which(carto$MASK_DWI_t0==1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

# display
resGR0 <- MRIaggr:::GRalgo(contrast=carto$TTP_t0, W=W, 
                           seed=indexN,sigma_max=1,range=c(7,+Inf), 
                           breaks=seq(-5,40,1),step=1,operator="sd",
                           iter_max=10,history_sigma=FALSE,history_step=TRUE,history_front=TRUE)

multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=3,hemisphere="lesion",legend=FALSE,
             breaks=seq(0,10,0.1),as.logical=TRUE,ylim=c(25,45), xlim=c(75,95),cex=2,            
             index1=list(coords=coords[resGR0$GR,],pch=20,cex=1),
             index2=list(coords=coords[indexN,],pch=20,cex=1)
)

multiplot(coords[resGR0$GR,],
             as.numeric(unlist(lapply(strsplit(resGR0$history_front,split=".",fixed=TRUE),"[",1)))
)


#### T2 
carto <- selectContrast(MRIaggr.Pat1_red,param=c("TTP_t0","MASK_DWI_t0","T2_FLAIR_t2"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,hemisphere="lesion")
W <- calcW(coords,range=sqrt(2))

indexN <- which(carto$MASK_DWI_t0==1)

multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",hemisphere="lesion")

system.time(
resGR0 <- MRIaggr:::GRalgo(contrast=carto$T2_FLAIR_t2, W=W, 
                           seed=indexN,sigma_max=500,range=c(300,+Inf),
                           breaks=seq(min(carto$T2_FLAIR_t2)-20,20+max(carto$T2_FLAIR_t2),by=20),
                           step=20,operator="sd",
                           iter_max=50,history_sigma=FALSE,history_step=TRUE,history_front=TRUE)
)
multiplot(coords,carto$T2_FLAIR_t2,
             index1=coords[resGR0$GR,])

index_tree <- 7
multiplot(coords[resGR0$GR,],             
             as.numeric(unlist(lapply(strsplit(resGR0$history_front,split=".",fixed=TRUE),"[",index_tree)))
)

tree <- strsplit(resGR0$history_front,split=".",fixed=TRUE)
stepTree <- sapply(1:resGR0$iter,function(x){as.numeric(unlist(lapply(tree,"[",x)))})
head(stepTree)

apply(stepTree,2,table)
table(resGR0$history_step)



#### 2- init ####
  
#### >a initGR - private #### 
# ?initGR
# MRIaggr:::initGR
# function(contrast,W,seed,range,range.seed,
# coords,sd.iter,breaks,scale,trace,method="initGR")
 
#### test
data(MRIaggr.Pat1_red, package="MRIaggr")

carto <- selectContrast(MRIaggr.Pat1_red,num=3,param=c("TTP_t0","MASK_DWI_t0"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,num=3,hemisphere="lesion")
W <- calcW(coords,range=sqrt(2))

indexN <- which(carto$MASK_DWI_t0==1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

res.init <- MRIaggr:::initGR(contrast=carto$TTP_t0,W=W,seed=seed,range=c(-Inf,Inf),range.seed=c(-Inf,Inf),
                             breaks=100,scale=FALSE,trace=TRUE,method="initGR")

#### example

