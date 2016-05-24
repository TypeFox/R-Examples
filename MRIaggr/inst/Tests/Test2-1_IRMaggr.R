#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Fichier test 2.1 : MRIaggr
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require(MRIaggr)

options(error=function() traceback(2)) 
options(max.print=10000)

#### 0 chargement ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 1- Select ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### >a selectContrast ####
# ?selectContrast
# findMethods("selectContrast",classes="MRIaggr")$MRIaggr
# .local <- function (object, param = NULL, num = NULL, format = "any", 
#  slice_var = "k", coords = FALSE, hemisphere = "both", 
# norm_mu = FALSE, norm_sigma = FALSE, na.rm = FALSE, subset = NULL) 

#### test
data("MRIaggr.Pat1_red", package="MRIaggr")
selectParameter(MRIaggr.Pat1_red)
sum(is.na(selectContrast(MRIaggr.Pat1_red,param="TTP_reperf")))
# multiplot(MRIaggr.Pat1_red,param="TTP_reperf")

length(selectContrast(MRIaggr.Pat1_red,param="TTP_reperf",na.rm=TRUE))

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select all parameters and all observations
carto <- selectContrast(MRIaggr.Pat1_red)
dim(carto)
head(carto)

## select a subset of parameters
carto <- selectContrast(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"))
dim(carto)
head(carto)

## select a subset of parameters on slices 1 to 3
carto <- selectContrast(MRIaggr.Pat1_red,num=1:3,param=c("DWI_t0","T2_FLAIR_t2"))
dim(carto)
head(carto)

## select a subset of parameters on slices 1 to 3 and normalized the center 
## the values using the contralateral
carto <- selectContrast(MRIaggr.Pat1_red,num=1:3,param=c("DWI_t0","T2_FLAIR_t2"),
                     norm_mu="contralateral")
dim(carto)
head(carto)

## select only observations which are lesioned at admission (i.e. MASK_DWI_t0=TRUE)
carto <- selectContrast(MRIaggr.Pat1_red,subset="MASK_DWI_t0",
                     param=c("DWI_t0","T2_FLAIR_t2","MASK_DWI_t0"))
dim(carto)
head(carto)

## select only observations which are lesioned at admission (i.e. MASK_DWI_t0=TRUE) with coordinates
carto <- selectContrast(MRIaggr.Pat1_red,subset="MASK_DWI_t0",
                     param=c("DWI_t0","T2_FLAIR_t2","MASK_DWI_t0"),coords=TRUE)
dim(carto)
head(carto)

## select only observations for which i=55
carto <- selectContrast(MRIaggr.Pat1_red,slice_var="i",num=55,coords=TRUE)
dim(carto)
head(carto)

#### >b selectClinic ####
# ?selectClinic
# findMethods("selectClinic",classes="MRIaggr")$MRIaggr
# .local <- function (object, param = NULL) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select all clinical data
res <- selectClinic(MRIaggr.Pat1_red)

## select only the gender
res <- selectClinic(MRIaggr.Pat1_red,param="Gender")

#### >c selectCoords ####
# ?selectCoords
# findMethods("selectCoords",classes="MRIaggr")$MRIaggr
# .local <- function (object, coords = c("i", "j", "k"), spatial_res = rep(1,3), num = NULL, 
#                     hemisphere = "both", subset = NULL, slice_var = "k",format = "data.frame") 
  
#### test
coords <- selectCoords(MRIaggr.Pat1_red,slice_var="j",num=16:18)
coords <- selectCoords(MRIaggr.Pat1_red,slice_var="j",num=16:18,subset="MASK_DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red,slice_var="j",num=50:70,subset="MASK_DWI_t0")

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select all coordinates for all observations
coords <- selectCoords(MRIaggr.Pat1_red)
dim(coords)
head(coords)

## select coordinate i for slices 4 and 6
coords <- selectCoords(MRIaggr.Pat1_red,coords="i",num=c(1,3))
dim(coords)
head(coords)

## select coordinate i for observations in the hemishere containing the lesion
coords <- selectCoords(MRIaggr.Pat1_red,hemisphere="lesion",num=c(1,3))
dim(coords)
head(coords)

## select coordinate i for observations in the right hemisphere
coords <- selectCoords(MRIaggr.Pat1_red,hemisphere="right",num=c(1,3))
dim(coords)
head(coords)

## select all coordinates and rescale them
coords <- selectCoords(MRIaggr.Pat1_red,spatial_res=c(1.875,1.875,6))
dim(coords)
head(coords)

## select coordinate i and j and return a matrix
coords <- selectCoords(MRIaggr.Pat1_red,format="matrix")
is(coords)
head(coords)


#### >d selectDefault_value  ####
# ?selectDefault_value
# findMethods("selectDefault_value",classes="MRIaggr")$MRIaggr
# .local <- function (object, param = NULL, as.numeric = FALSE)

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select default values for all parameters
res <- selectDefault_value(MRIaggr.Pat1_red)

## select default values for T2_FLAIR_t2 and DWI_t0
res <- selectDefault_value(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","DWI_t0"))

## select default values for T2_FLAIR_t2 and DWI_t0 and convert them in numeric
res <- selectDefault_value(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","DWI_t0"),as.numeric=TRUE)

#### >e selectDescStats ####
# ?selectDescStats
# findMethods("selectDescStats",classes="MRIaggr")$MRIaggr
# .local <- function (object, name = NULL, subset_W = NULL, hemisphere = "both", num = NULL, slice_var = "k") 

#### test
  
#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select all elements in the slot @ls_descStats
ls_descStats <- selectDescStats(MRIaggr.Pat1_red)
names(ls_descStats)

## get the name of all elements present in the slot @ls_descStats
selectParameter(MRIaggr.Pat1_red,type="ls_descStats")

## select a specific element 
res <- selectDescStats(MRIaggr.Pat1_red,name="GroupsLesion")

## compute and affect a neighborhood matrix
calcW(MRIaggr.Pat1_red,range=3,update.object=TRUE,overwrite=TRUE,
      spatial_res=c(1.875,1.875,6))

## select the neighborhood matrix for a subset of observations
res <- selectDescStats(MRIaggr.Pat1_red,name="W_euclidean",hemisphere="lesion",num=3)
dim(res)
table(rowSums(res>0))

res <- selectDescStats(MRIaggr.Pat1_red,name="W_euclidean",subset_W = 1:10)
res

#### >f selectVoxelDim  ####
# ?selectVoxelDim
# findMethods("selectVoxelDim",classes="MRIaggr")$MRIaggr
# .local <- function (object) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## selection
selectVoxelDim(MRIaggr.Pat1_red)

#### >g selectHistory  ####
# ?selectHistory
# findMethods("selectHistory",classes="MRIaggr")$MRIaggr
# .local <- function (object) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## selection
selectHistory(MRIaggr.Pat1_red)

#### >h selectIdentifier ####
# ?selectIdentifier
# findMethods("selectIdentifier",classes="MRIaggr")$MRIaggr
# .local <- function (object) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## selection
selectIdentifier(MRIaggr.Pat1_red)

#### >i selectN ####
# ?selectN
# findMethods("selectN",classes="MRIaggr")$MRIaggr
# .local <- function (object, num = NULL, hemisphere = "both", subset = NULL) 

#### test

#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## total number of observations
res <- selectN(MRIaggr.Pat1_red)

## number of observations for the hemisphere that contains the lesion
res <- selectN(MRIaggr.Pat1_red,hemisphere="lesion")

## number of observations in the first 1000 observations that are in the  hemisphere that contains the lesion
res <- selectN(MRIaggr.Pat1_red,subset=1:1000,hemisphere="lesion")

#### >j selectNormalization ####
# ?selectNormalization
# findMethods("selectNormalization",classes="MRIaggr")$MRIaggr
# .local <- function (object, type = NULL, mu = TRUE, sigma = TRUE, 
#                     hemisphere = "both", num = NULL, param = NULL) 
  
#### test
res <- selectNormalization(MRIaggr.Pat1_red,type="global",mu=TRUE,sigma=FALSE,hemisphere="left")
res <- selectNormalization(MRIaggr.Pat1_red,type="slice",mu=TRUE,sigma=FALSE,num=1,param="DWI_t0")

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select all normalization values
res <- selectNormalization(MRIaggr.Pat1_red)
names(res)

## select specific normalization normalization values 
# computed on the whole brain
res <- selectNormalization(MRIaggr.Pat1_red,type="global",mu=TRUE,sigma=FALSE,hemisphere="both")
# idem but only for DWI_t0
res <- selectNormalization(MRIaggr.Pat1_red,type="global",mu=TRUE,sigma=FALSE,param="DWI_t0") 

# computed by slice
res <- selectNormalization(MRIaggr.Pat1_red,type="slice",mu=TRUE,sigma=FALSE,hemisphere="both") 
# idem for slice 1
res <- selectNormalization(MRIaggr.Pat1_red,type="slice",mu=TRUE,sigma=FALSE,num=1) 

# computed on 3 consecutive slices
res <- selectNormalization(MRIaggr.Pat1_red,type="3slices",mu=FALSE,sigma=TRUE, 
                           hemisphere="left",num=2,param="T2_FLAIR_t2")

#### >k selectParameter ####
# ?selectParameter
# findMethods("selectParameter",classes="MRIaggr")$MRIaggr
# .local <- function (object, type = "contrast", mask = TRUE) 

#### test

#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## extract all imaging parameters
res <- selectParameter(MRIaggr.Pat1_red)

## extract the slot @hemiLesion
res <- selectParameter(MRIaggr.Pat1_red,type="hemiLesion")

## extract the names of the parameters in the slot @ls_descStats
res <- selectParameter(MRIaggr.Pat1_red,type="ls_descStats")

#### >l selectTable ####
# ?selectTable
# findMethods("selectTable",classes="MRIaggr")$MRIaggr
# .local <- function (object, type) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## selection
res <- selectTable(MRIaggr.Pat1_red,type="lesion")
res <- selectTable(MRIaggr.Pat1_red,type="reperfusion")
res <- selectTable(MRIaggr.Pat1_red,type="hypoperfusion")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 2- affect ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### >a affectContrast ####
# ?"affectContrast<-"
# findMethods("affectContrast<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, param = NULL, default_value = NULL, 
#                     overwrite = FALSE, trace = TRUE) 

#### test
affectContrast(MRIaggr.Pat1_red,param="noise",default_value=data.frame("noise"=0),overwrite=TRUE) <- rnorm(selectN(MRIaggr.Pat1_red))

affectContrast(MRIaggr.Pat1_red,param="mask",default_value=data.frame("mask"=0),overwrite=TRUE) <- rnorm(selectN(MRIaggr.Pat1_red))>10

#### example 
## load NIFTI files and convert them to MRIaggr
path <- system.file(file.path("nifti"),package = "MRIaggr")
ls.array <- list()
ls.array[[1]] <- readMRI(file=file.path(path,"DWI_t0"),format="nifti")
ls.array[[2]] <- readMRI(file=file.path(path,"MASK_DWI_t0"),format="nifti")
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=c("DWI_t0","MASK_DWI_t0"))

## affect a new contrast parameters
affectContrast(MRIaggr.Pat1,param="noise",overwrite=TRUE) <- rnorm(selectN(MRIaggr.Pat1))

## perform operations on a contrast parameters and store the results
myCarto <- selectContrast(MRIaggr.Pat1,param="DWI_t0")
myCarto <- myCarto*2+1
affectContrast(MRIaggr.Pat1,param="myCarto",overwrite=TRUE) <- myCarto

## import a contrast parameters in an already existing MRIaggr object
nifti.MTT_t0 <- readMRI(file=file.path(path,"MTT_t0"),format="nifti")
df.MTT_t0 <- array2df(nifti.MTT_t0,name_newparam="MTT_t0")$MTT_t0
affectContrast(MRIaggr.Pat1,param="MTT_t0",overwrite=TRUE) <- df.MTT_t0

## some calc methods automatically save results in the @data slot
calcFilter(MRIaggr.Pat1,param="MTT_t0",filter="2D_G3",
           update.object=TRUE,overwrite=TRUE)
res <- selectContrast(MRIaggr.Pat1,param="MTT_t0_2D_G3")

#### >b affectClinic  ####
# ?"affectClinic<-"
# findMethods("affectClinic<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, add = FALSE, overwrite = FALSE, trace = TRUE) 

#### test

# affect empty clinical data
affectClinic(MRIaggr.Pat1_red) <- data.frame()

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## affect clinical data
affectClinic(MRIaggr.Pat1_red) <- data.frame(Age=32,Gender="Male",NIHSS_H0="5")
selectClinic(MRIaggr.Pat1_red,param="Age")

## add a new parameter
affectClinic(MRIaggr.Pat1_red,add=TRUE,overwrite=TRUE) <- data.frame(City="Lyon")
selectClinic(MRIaggr.Pat1_red)

#### >c affectDescStats  ####
# ?"affectDescStats<-"
# findMethods("affectDescStats<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, name, overwrite = FALSE, trace = TRUE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## affect a vector
affectDescStats(MRIaggr.Pat1_red,name="spatial_res") <- c(1.875,1.875,6)

## select the corresponding element
selectDescStats(MRIaggr.Pat1_red,"spatial_res")

## some calc methods automatically save results in the ls_descStats slot
# find spatial groups 
calcGroupsMask(MRIaggr.Pat1_red,mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               bandwidth=6,
               spatial_res=selectDescStats(MRIaggr.Pat1_red,"spatial_res"),
               update.object=TRUE,overwrite=TRUE)

# extract spatial groups
selectDescStats(MRIaggr.Pat1_red,"GroupsLesion")


#### >d affectHemisphere ####
# ?"affectHemisphere<-"
# findMethods("affectHemisphere<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, overwrite = FALSE, trace = TRUE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## manual affectation
resHemi <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,i_test=2,j_test=2,angle_test=2)
affectHemisphere(MRIaggr.Pat1_red,overwrite=TRUE) <- list(midplane=resHemi$midplane,
                                                          data=resHemi$data)

## automatic affectation
resHemi <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,i_test=2,j_test=2,angle_test=2,
                          update.object=TRUE,overwrite=TRUE)

## display
index1 <- data.frame(selectParameter(MRIaggr.Pat1_red,type="midplane"),15)
names(index1) <- c("i","j","k")
multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,midplane=TRUE,window=FALSE,
             index1=list(coords=index1,pch=20,cex=2,col="purple")
             
)

#### >e affectNormalization ####
# ?"affectNormalization<-"
# findMethods("affectNormalization<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, overwrite = FALSE, trace = TRUE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## parameters to normalize
param <- c("DWI_t0","T2_FLAIR_t2","T2_GRE_t0","TTP_t0")

## manual affectation
resNormalization <- calcNormalization(MRIaggr.Pat1_red,param=param)
affectNormalization(MRIaggr.Pat1_red,overwrite=TRUE) <- resNormalization

## automatic affectation
resNormalization <- calcNormalization(MRIaggr.Pat1_red,param=param,
                                      update.object=TRUE,overwrite=TRUE)

#### >f affectTable ####
# ?"affectTable<-"
# findMethods("affectTable<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, type, overwrite = FALSE, 
#                     trace = TRUE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

#### 1- lesion 
## manual affectation
maskN <- c("MASK_DWI_t0","MASK_T2_FLAIR_t2")
resTable <- calcTableLesion(MRIaggr.Pat1_red,maskN=maskN,as.logical=TRUE)
affectTable(MRIaggr.Pat1_red,type="lesion",overwrite=TRUE) <- resTable

## automatic affectation
resTable <- calcTableLesion(MRIaggr.Pat1_red,maskN=maskN,
                            as.logical=TRUE,update.object=TRUE,overwrite=TRUE)

## display
selectTable(MRIaggr.Pat1_red,type="lesion")

#### 2- hypoperfusion and reperfusion 
## manual affectation
resTable <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP","MTT"),time=c("t0","t1"))
affectTable(MRIaggr.Pat1_red,type="hypoperfusion",overwrite=TRUE) <- resTable$volume_hypo
affectTable(MRIaggr.Pat1_red,type="reperfusion",overwrite=TRUE) <- resTable$volume_reperf

## automatic affectation
resTable <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP","MTT"),time=c("t0","t1"),
                                update.object=TRUE,overwrite=TRUE)

## display
selectTable(MRIaggr.Pat1_red,type="reperfusion")

#### >g supprContrast  ####
# ?"supprContrast<-"
# findMethods("supprContrast<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, trace = TRUE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## available contrast parameters
selectParameter(MRIaggr.Pat1_red)

## delete two contrast parameters
supprContrast(MRIaggr.Pat1_red) <- c("MTT_t0","MASK_DWI_t0")

## remaining contrast parameters
selectParameter(MRIaggr.Pat1_red)

#### >h supprDescStats  ####
# ?"supprDescStats<-"
# findMethods("supprDescStats<-",classes="MRIaggr")$MRIaggr
# .local <- function (object, value, trace = TRUE) 

#### test
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## existing elements in @ls_descStats
selectParameter(MRIaggr.Pat1_red,type="ls_descStats")

## delete one element in @ls_descStats
supprDescStats(MRIaggr.Pat1_red) <- "index_sauve"

## remaining elements in @ls_descStats
selectParameter(MRIaggr.Pat1_red,"ls_descStats")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 3- calc ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### >a calcBrainMask ####
# ?calcBrainMask
# findMethods("calcBrainMask",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, type = "kmeans", th.nb_breaks = 100, 
#                     th.nb_optima = 5, th.select_optima = 1, th.upper = TRUE, 
#                     kmeans.n_groups = 2, kmeans.Neighborhood = 3, skull.param = NULL, 
#                     skull.n_groups = 3, plot = TRUE, window = FALSE, 
#                     filename = "auto", width = 1000, height = 700, path = NULL, 
#                     unit = "px", res = NA, trace = TRUE, update.object = FALSE, 
#                     overwrite = FALSE) 

#### test

#### example 
## load NIFTI files and convert them to MRIaggr
path <- system.file(file.path("nifti"),package = "MRIaggr")
ls.array <- list(readMRI(file=file.path(path,"T1_t0"),format="nifti"),
                 readMRI(file=file.path(path,"T2_GRE_t0"),format="nifti"))
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=c("T1_t0","T2_GRE_t0"))

#### 1- thresholding approach 
res <- calcBrainMask(MRIaggr.Pat1,param="T2_GRE_t0",type="threshold",
                     th.select_optima=2)

breaks <- res$analysis[,"threshold"]
res <- calcBrainMask(MRIaggr.Pat1,param="T2_GRE_t0",type="threshold",
                     th.breaks=breaks[breaks>50],th.select_optima=1,
                     overwrite=TRUE,update.object=TRUE)

## display
multiplot(MRIaggr.Pat1,param="mask")

multiplot(MRIaggr.Pat1,param="T2_GRE_t0",index1="mask")

## other parameter 
res <- calcBrainMask(MRIaggr.Pat1,param="T1_t0",type="threshold",
                     th.breaks=200)

res <- calcBrainMask(MRIaggr.Pat1,param="T1_t0",type="threshold",
                     th.breaks=seq(0,400,length.out=50),th.select_optima=2,
                     overwrite=TRUE,update.object=TRUE)

multiplot(MRIaggr.Pat1,param="mask")


#### 2- k-means approach 
res <- calcBrainMask(MRIaggr.Pat1,param="T2_GRE_t0",type="kmeans",
                     kmeans.n_groups=2:4,
                     update.object=TRUE,overwrite=TRUE)

## display
multiplot(MRIaggr.Pat1,param="T2_GRE_t0",index1="mask")
multiplot(MRIaggr.Pat1,param="mask")


#### >b calcContralateral ####
# ?calcContralateral
# findMethods("calcContralateral",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, num = NULL, type = "mean", 
#                     param.ref = NULL, distband = 1, lambda = 1, trace = TRUE, 
#                     update.object = FALSE, overwrite = FALSE) 

#### test
res <- calcContralateral(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),num=NULL,type="mean",trace=TRUE)
multiplot(res$data[,c("i","j","k")],
             res$data[,"DWI_t0_contro",drop=FALSE])
res <- calcContralateral(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),num=NULL,type="median",trace=TRUE)
multiplot(res$data[,c("i","j","k")],
             res$data[,"DWI_t0_contro",drop=FALSE])
res <- calcContralateral(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),num=NULL,type="1NN_penalised",param.ref="T1_t0",
                         distband=1,lambda=1,trace=TRUE)
multiplot(res$data[,c("i","j","k")],
             res$data[,"DWI_t0_contro",drop=FALSE])

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## associate each observation to its contralateral correspondant 
## according T1 parameter and compute the normalized parameters
res <- calcContralateral(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),num=NULL,
                         type="mean",param.ref="T1_t0", distband=1,lambda=1,trace=TRUE)

## display
multiplot(res$data[,c("i","j","k"),drop=FALSE],
             contrast=res$data$DWI_t0_contro
)

multiplot(res$data[,c("i","j","k"),drop=FALSE],
             contrast=res$data$DWI_t0_contro,
             index1=res$data[res$index_plot$index.plot_lesionR,c("i","j","k"),drop=FALSE],
             index2=res$data[res$index_plot$index.plot_lesionL,c("i","j","k"),drop=FALSE]
)

#### >c calcDistMask ####
# ?calcDistMask
# findMethods("calcDistMask",classes="MRIaggr")$MRIaggr
# .local <- function (object, mask, name_newparam = NULL, spatial_res = c(1,1, 1), 
# as.logical = FALSE, filter = "3D_N10", trace = TRUE, 
# update.object = FALSE, overwrite = FALSE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute neighborhood matrix 
MASK_DWI_t0 <- selectContrast(MRIaggr.Pat1_red,param="MASK_DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)
W <- calcW(object=as.data.frame(coords[MASK_DWI_t0==1,]),range=sqrt(2),row.norm=TRUE)

## compute spatial groups
res.Groups <- calcGroupsW(W)
res.Groups$group_size

## display
multiplot(coords[MASK_DWI_t0==1,],contrast=res.Groups$group,
             legend=FALSE,xlim=c(30,100),ylim=c(20,100),cex=0.5,
             index1=list(coords=coords,outline=TRUE))

#### >d calcDistTissues ####
# ?calcDistTissues
# findMethods("calcDistTissues",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, class, num = NULL, hemisphere = "both") 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute the distribution of DWI and T2 FLAIR for the CSF, WM, GM and lesion observations
res <- calcDistTissues(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),
                       class=c("CSF","WM","GM","MASK_DWI_t0")
)

#### >e calcFilter ####
# ?calcFilter
# findMethods("calcFilter",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, filter, w_contrast = FALSE, 
#                     na.rm = FALSE, name_newparam = NULL, trace = TRUE, update.object = FALSE, 
#                     overwrite = FALSE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute and affect filtered parameter to the MRIaggr object
# gaussian filter
calcFilter(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","DWI_t0","TTP_t0"),
           filter="2D_G3",w_contrast=FALSE,na.rm=FALSE,update.object=TRUE,overwrite=TRUE)
selectParameter(MRIaggr.Pat1_red)

# median filter
calcFilter(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","DWI_t0","TTP_t0"),
           filter="2D_M3",na.rm=FALSE,update.object=TRUE,overwrite=TRUE)

## display
par(mfrow=c(2,2))
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
             num=1,window=NULL,breaks=c(-100,seq(0,450),601))
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2_2D_G3",
             num=1,legend=FALSE,window=NULL,breaks=c(-100,seq(0,450),601))
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2_2D_M3",
             num=1,legend=FALSE,window=NULL,breaks=c(-100,seq(0,450),601))

## see the results of the different filters
# G : Gaussian filter
resG <- calcFilter(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
                   filter="2D_G3")

# M : median filter
resM <- calcFilter(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
                   filter="2D_M3")

# S : Sobel filter
resS <- calcFilter(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
                   filter="2D_Sx")

# I
resI.T <- calcFilter(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",
                     filter="2D_I3",norm=TRUE)
resI.F <- calcFilter(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",
                     filter="2D_I3",norm=FALSE)

# N
resN.T <- calcFilter(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",
                     filter="3D_N10",norm=TRUE)
resN.F <- calcFilter(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",
                     filter="3D_N10",norm=FALSE)

## display
par(mfrow=c(2,2),mar=rep(2,4),mgp=c(2,0.75,0))
breaks <- seq(-50,500,1)
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=3,
             breaks=breaks,window=NULL,legend=FALSE,
             main="no filtering",num.main=FALSE,main.legend="")
multiplot(resS$res[,c("i","j","k")],contrast=resS$res[,4],
             num=3,window=NULL,legend=FALSE,
             palette="cm.colors",breaks=seq(-1000,1000),
             main="sobelX filtering",num.main=FALSE)
multiplot(resG$res[,c("i","j","k")],contrast=resG$res[,4],
             num=3,window=NULL,legend=FALSE,breaks=breaks,
             main="gaussian filtering",num.main=FALSE)
multiplot(resM$res[,c("i","j","k")],contrast=resM$res[,4],
             num=3,window=NULL,legend=FALSE,breaks=breaks,
             main="median filtering",num.main=FALSE)

layout(matrix(1:6,nrow=3,ncol=2,byrow=TRUE),widths=c(2,1),heights=rep(1,3))
par(mar=rep(2,4),mgp=c(2,0.75,0))
multiplot(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",num=3,
             window=NULL,legend=TRUE,main="raw",num.main=FALSE,main.legend="")
multiplot(resI.T$res[,c("i","j","k")],contrast=resI.T$res[,4],
             num=3,window=NULL,legend=TRUE,
             main="Influence filtering - norm=TRUE",num.main=FALSE)
multiplot(resI.F$res[,c("i","j","k")],contrast=resI.F$res[,4],
             num=3,window=NULL,legend=TRUE,
             main="Influence filtering - norm=FALSE",num.main=FALSE)


layout(matrix(1:6,nrow=3,ncol=2,byrow=TRUE),widths=c(2,1),heights=rep(1,3))
par(mar=rep(2,4),mgp=c(2,0.75,0))
multiplot(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",num=3,
             window=NULL,legend=TRUE,main="raw",num.main=FALSE,main.legend="")
multiplot(resN.T$res[,c("i","j","k")],contrast=resN.T$res[,4],
             num=3,window=NULL,legend=TRUE,
             main="Neighborhood3D filtering - norm=TRUE",num.main=FALSE)
multiplot(resN.F$res[,c("i","j","k")],contrast=resN.F$res[,4],
             num=3,window=NULL,legend=TRUE,
             main="Neighborhood3D filtering - norm=FALSE",num.main=FALSE)

#### >f calcGroupsMask ####
# ?calcGroupsMask
# findMethods("calcGroupsMask",classes="MRIaggr")$MRIaggr
# .local <- function (object, mask, bandwidth, spatial_res = c(1, 1, 1), 
#                     as.logical = FALSE, W = "ifany", trace = TRUE, 
#                     update.object = FALSE, overwrite = TRUE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

# compute spatial groups
calcGroupsMask(MRIaggr.Pat1_red,mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               W.range=6,W.spatial_res=c(1.875,1.875,6),
               update.object=TRUE,overwrite=TRUE)

# extract spatial groups
selectDescStats(MRIaggr.Pat1_red,"GroupsLesion")

#### >g calcHemisphere ####
# ?calcHemisphere
# findMethods("calcHemisphere",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, num = NULL, p = 1, subset = NULL, 
#                     penalty = "symmetry", mask = NULL, as.logical = FALSE, 
#                     i_test = 5, j_test = 5, angle_test = 5, unit_angle = "radian", 
#                     n.points = 100, plot = TRUE, window = FALSE, filename = "auto", 
#                     width = 1000, height = 700, path = NULL, unit = "px", 
#                     res = NA, trace = TRUE, update.object = FALSE, overwrite = FALSE) 

#### test
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,p=1,penalty="asymmetry",
                      i_test=2,j_test=2,angle_test=2,
                      update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
             midplane=TRUE)

res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,p=1,penalty="symmetry",
                      i_test=2,j_test=2,angle_test=2,n.points=50)

res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,p=1,penalty="symmetry",
                      i_test=2,j_test=2,angle_test=2,n.points=50,plot=FALSE)


#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,
                      i_test=2,j_test=2,angle_test=2,
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)

## compute the mid-sagittal plan
multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0",num=3,legend=FALSE,
             midplane=TRUE,main="original coordinates - slice ")

## display
multiplot(selectContrast(MRIaggr.Pat1_red,param=c("i_hemisphere","j_hemisphere","k")),
             contrast=selectContrast(MRIaggr.Pat1_red,param="T2_GRE_t0"),num=3,
             index1=cbind(0,seq(-50,50),3),main="new coordinates - slice ",legend=FALSE)


## compute the mid-sagittal plan and mark lesion
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,
                      mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), i_test=2,j_test=2,angle_test=2,
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)

#### >h calcNormalization ####
# ?calcNormalization
# findMethods("calcNormalization",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, mu_type = "mean", sigma_type = "sd", 
#                     rm.CSF = FALSE, rm.WM = FALSE, rm.GM = FALSE, 
#                     trace = TRUE, update.object = FALSE, overwrite = FALSE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute normalization values
res <- calcNormalization(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),
                         update.object=TRUE,overwrite=TRUE)

## display
par(mfrow=c(2,4),mar=rep(1.5,4),mgp=c(2,0.5,0))
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=1:3,
             legend=TRUE,window=NULL,main="raw - slice ")
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=1:3,
             norm_mu="contralateral",norm_sigma="contralateral",
             legend=TRUE,window=NULL,main="normalized - slice ")

## extract normalization
selectNormalization(MRIaggr.Pat1_red,type="global",mu=TRUE,sigma=FALSE)

#### >i calcRegionalContrast ####
# ?calcRegionalContrast
# findMethods("calcRegionalContrast",classes="MRIaggr")$MRIaggr


#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute regional values
res  <- calcRegionalContrast(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","T1_t0"),
                        W.spatial_res=c(1.875,1.875,1.875),W.range=6,bandwidth=1.875,
                        update.object=TRUE,overwrite=TRUE)

## display
par(mfrow=c(2,4),mar=rep(1.5,4),mgp=c(2,0.5,0))
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=1:3,
             window=NULL,main="raw - slice ")
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2_regional",num=1:3,
             window=NULL,main="regional - slice ")

#### >j calcROCthreshold ####
# ?calcROCthreshold
# findMethods("calcROCthreshold",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, mask, as.logical = FALSE, 
#                     digit = 2, plot = "ROC", digit.plot = 3, window = FALSE, 
#                     filename = "auto", width = 1000, height = 700, path = NULL, 
#                     unit = "px", res = NA, trace = TRUE, update.object = FALSE, 
#                     overwrite = FALSE) 
#### test
res <- calcROCthreshold(MRIaggr.Pat1_red,mask=c("MASK_DWI_t0"),param=c("DWI_t0"),
                        as.logical=TRUE)

res <- calcROCthreshold(MRIaggr.Pat1_red,mask=c("MASK_DWI_t0"),param=c("DWI_t0"),
                        as.logical=TRUE,plot="boxplot_Youden")

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## ROC analysis
res <- calcROCthreshold(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),
                        mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),as.logical=TRUE)

res <- calcROCthreshold(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),
                        mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),as.logical=TRUE,
                        plot="boxplot_Youden")

#### >k calcSmoothMask ####
# ?calcSmoothMask
# findMethods("calcSmoothMask",classes="MRIaggr")$MRIaggr
# .local <- function (object, mask = "mask", as.logical = FALSE, 
#                     size_2Dgroup = 50, Neighborhood_2D = "3D_N8", rm.2Dhole = FALSE, 
#                     size_3Dgroup = "unique", Neighborhood_3D = "3D_N10", 
#                     rm.3Dhole = TRUE, erosion_th = 0.75, Vmask_min = 0.25, 
#                     Vbackground_max = 0.75, Neighborhood_V = "3D_n10", trace = TRUE, 
#                     update.object = FALSE, overwrite = FALSE) 

#### test
res <- calcSmoothMask(MRIaggr.Pat1,rm.2Dhole = TRUE, 
                      update.object=TRUE,overwrite=TRUE)

#### example 
## load data and build MRIaggr
path <- system.file(file.path("nifti"),package = "MRIaggr")
ls.array <- list(readMRI(file=file.path(path,"T2_GRE_t0"),format="nifti"))
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param="T2_GRE_t0")

## create the cerebral mask
res <- calcBrainMask(MRIaggr.Pat1,param="T2_GRE_t0",type="kmeans",
                     kmeans.n_groups=2:4,
                     update.object=TRUE,overwrite=TRUE)

## smooth the cerebral mask
res <- calcSmoothMask(MRIaggr.Pat1,update.object=TRUE,overwrite=TRUE)

## display
multiplot(MRIaggr.Pat1,param="mask",legend=FALSE)

#### >l calcTableHypoReperf ####
# ?calcTableHypoReperf
# findMethods("calcTableHypoReperf",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, time, threshold = 1:10, 
#                     sep = "_", rm.CSF = FALSE, mask = NULL, as.logical = FALSE, 
#                     hemisphere = "both", norm_mu = FALSE, norm_sigma = FALSE, 
#                     subset = NULL, trace = TRUE, param.update = "reperf", 
#                     update.object = FALSE, overwrite = FALSE) 
  
#### test
carto_TTP_t0 <- selectContrast(MRIaggr.Pat1_red,param="TTP_t0")
carto_TTP_t1 <- selectContrast(MRIaggr.Pat1_red,param="TTP_t1")

res <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP","MTT"),time=c("t0","t1"),
                           mask="MASK_DWI_t0",as.logical=TRUE,
                           update.object=TRUE,overwrite=TRUE)

# hypo
for(test.threshold in c(1,4,10)){
  tempo1 <- sum( (carto_TTP_t0>=test.threshold) )-res$volume_hypo[as.character(test.threshold),"Vhypo.TTP_t0"]
  tempo2 <- sum( (carto_TTP_t1>=test.threshold) )-res$volume_hypo[as.character(test.threshold),"Vhypo.TTP_t1"]
  cat(test.threshold," : ",tempo1," ",tempo2,"\n")
}

# mismatch
testN <- (selectContrast(MRIaggr.Pat1_red,param="MASK_DWI_t0")==0)
for(test.threshold in c(1,4,10)){
  tempo1 <- sum( (carto_TTP_t0>=test.threshold)*testN )-res$volume_hypo[as.character(test.threshold),"Vmismatch.TTP"]
  tempo2 <- sum( (carto_TTP_t0>=test.threshold)*testN )/sum(testN==FALSE)-res$volume_hypo[as.character(test.threshold),"PCmismatch.TTP"]
  cat(test.threshold," : ",tempo1," ",tempo2,"\n")
}

# reperf
test.threshold <- 1
for(test.threshold in c(1,4,10)){
  tempo1 <- sum((carto_TTP_t0>=test.threshold)*(carto_TTP_t1<test.threshold))-res$volume_reperf[as.character(test.threshold),"Vreperf.TTP"]
  tempo2 <- sum((carto_TTP_t0>=test.threshold)*(carto_TTP_t1<test.threshold))/sum(carto_TTP_t0>=test.threshold)-res$volume_reperf[as.character(test.threshold),"PCreperf.TTP"]
  tempo3 <- sum((carto_TTP_t0<test.threshold)*(carto_TTP_t1>=test.threshold))-res$volume_reperf[as.character(test.threshold),"Vdeperf.TTP"]
  tempo4 <- sum((carto_TTP_t0<test.threshold)*(carto_TTP_t1>=test.threshold))/sum(carto_TTP_t0>=test.threshold)-res$volume_reperf[as.character(test.threshold),"PCdeperf.TTP"]
  
  
  cat(test.threshold," : ",tempo1," ",tempo2," ",tempo3," ",tempo4,"\n")
}

# Wreperf
carto_TTP_t0.norm <- carto_TTP_t0
carto_TTP_t0.norm[carto_TTP_t0.norm<0] <- 0
carto_TTP_t0.norm[carto_TTP_t0.norm>10] <- 10

carto_TTP_t1.norm <- carto_TTP_t1
carto_TTP_t1.norm[carto_TTP_t1.norm<0] <- 0
carto_TTP_t1.norm[carto_TTP_t1.norm>10] <- 10

# weight
shift_TTP.norm <- carto_TTP_t1.norm-carto_TTP_t0.norm
index.reperf <- which((shift_TTP.norm<0)*(carto_TTP_t0.norm>0)==1)
index.deperf <- which((shift_TTP.norm>0)*(carto_TTP_t0.norm<10)==1)
w_reperf <- rep(0,length(carto_TTP_t0.norm))
w_reperf[index.reperf] <- -shift_TTP.norm[index.reperf]/carto_TTP_t0.norm[index.reperf]
# boxplot(w_reperf[w_reperf>0])
w_deperf <- rep(0,length(carto_TTP_t0.norm))
w_deperf[index.deperf] <- shift_TTP.norm[index.deperf]/(10-carto_TTP_t0.norm[index.deperf])
# boxplot(w_deperf[w_deperf>0])

for(test.threshold in c(1,4,10)){
  tempo1 <- sum(w_reperf*(carto_TTP_t0>=test.threshold)*(carto_TTP_t1<test.threshold))-res$volume_reperf[as.character(test.threshold),"VreperfW.TTP"]
  tempo2 <- sum(w_reperf*(carto_TTP_t0>=test.threshold)*(carto_TTP_t1<test.threshold))/sum(carto_TTP_t0>=test.threshold)-res$volume_reperf[as.character(test.threshold),"PCreperfW.TTP"]
  tempo3 <- sum(w_deperf*(carto_TTP_t0<test.threshold)*(carto_TTP_t1>=test.threshold))-res$volume_reperf[as.character(test.threshold),"VdeperfW.TTP"]
  tempo4 <- sum(w_deperf*(carto_TTP_t0<test.threshold)*(carto_TTP_t1>=test.threshold))/sum(carto_TTP_t0>=test.threshold)-res$volume_reperf[as.character(test.threshold),"PCdeperfW.TTP"]
  
  cat(test.threshold," : ",tempo1," ",tempo2," ",tempo3," ",tempo4,"\n")
}

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

#### 1- directly
res <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP","MTT"),time=c("t0","t1"),
                           mask="MASK_DWI_t0",as.logical=TRUE,
                           update.object=TRUE,overwrite=TRUE)

carto_TTP_t0 <- selectContrast(MRIaggr.Pat1_red,param="TTP_t0")
carto_TTP_t1 <- selectContrast(MRIaggr.Pat1_red,param="TTP_t1")

## hypoperfusion
sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red,"hypoperfusion")["4","Vhypo.TTP_t0"]

## mismatch
testN <- (selectContrast(MRIaggr.Pat1_red,param="MASK_DWI_t0")==0)

sum( (carto_TTP_t0>=4)*testN )
selectTable(MRIaggr.Pat1_red,"hypoperfusion")["4","Vmismatch.TTP"]

sum( (carto_TTP_t0>=4)*testN )/sum( testN==FALSE )
selectTable(MRIaggr.Pat1_red,"hypoperfusion")["4","PCmismatch.TTP"]


## reperfusion
sum((carto_TTP_t0>=4)*(carto_TTP_t1<4))
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","Vreperf.TTP"]

sum((carto_TTP_t0>=4)*(carto_TTP_t1<4))/sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","PCreperf.TTP"]

## W reperfusion
carto_TTPth_t0 <- carto_TTP_t0
carto_TTPth_t0[carto_TTPth_t0>10] <- 10
carto_TTPth_t0[carto_TTPth_t0<0] <- 0

carto_TTPth_t1 <- carto_TTP_t1
carto_TTPth_t1[carto_TTP_t1>10] <- 10
carto_TTPth_t1[carto_TTP_t1<0] <- 0

weight <- (carto_TTPth_t0-carto_TTPth_t1)/carto_TTPth_t0
weight[ ((carto_TTPth_t0==0)+(carto_TTP_t0<4)+(carto_TTP_t1>=4)) > 0 ] <- 0

sum((carto_TTP_t0>=4)*(carto_TTP_t1<4)*weight)
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","VreperfW.TTP"]

sum((carto_TTP_t0>=4)*(carto_TTP_t1<4)*weight)/sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","PCreperfW.TTP"]


## deperfusion
sum((carto_TTP_t0<4)*(carto_TTP_t1>=4))
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","Vdeperf.TTP"]

sum((carto_TTP_t0<4)*(carto_TTP_t1>=4))/sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","PCdeperf.TTP"]

## shift
sum((carto_TTPth_t0-carto_TTPth_t1>=4))
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","Vshift_reperf.TTP"]

sum((carto_TTPth_t0-carto_TTPth_t1>=4))/sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red,"reperfusion")["4","PCshift_reperf.TTP"]

#### 2- via calcThresholdMRIaggr 
calcThresholdMRIaggr(MRIaggr.Pat1_red,param=c("TTP_t0","MTT_t0","TTP_t1","MTT_t1"),
                 threshold=1:10,name_newparam=c("TTP.GR_t0","MTT.GR_t0","TTP.GR_t1","MTT.GR_t1"),
                 rm.CSF=TRUE,hemisphere="lesion",
                 GRalgo=TRUE,seed=c("MASK_T2_FLAIR_t2","MASK_DWI_t0"),W=NULL,bandwidth=sqrt(2),
                 update.object=TRUE,overwrite=TRUE)

res <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP.GR","MTT.GR"),time=c("t0","t1"),
                           mask="MASK_DWI_t0",as.logical=TRUE,
                           update.object=TRUE,overwrite=TRUE)

## display
selectTable(MRIaggr.Pat1_red,"hypoperfusion")["4","Vhypo.TTP.GR_t0"]

par(mfrow=c(2,4),mar=rep(1.5,4),mgp=c(2,0.5,0))
multiplot(MRIaggr.Pat1_red,param="TTP_t0",num=1:3,
             palette=rainbow(10),window=NULL,main="raw - slice ",breaks=(0:10)-10^{-10})
multiplot(MRIaggr.Pat1_red,param="TTP.GR_t0",num=1:3,
             palette=rainbow(10),window=NULL,main="GR - slice ",breaks=(0:10)-10^{-10})


#### >m calcTableLesion ####
# ?calcTableLesion
# findMethods("calcTableLesion",classes="MRIaggr")$MRIaggr
# .local <- function (object, maskN, mask = NULL, as.logical = FALSE, 
#                     trace = TRUE, update.object = FALSE, overwrite = FALSE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute table
res <- calcTableLesion(MRIaggr.Pat1_red,maskN=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
                       as.logical=TRUE,update.object=TRUE,overwrite=TRUE)

## extract table
res <- selectTable(MRIaggr.Pat1_red,"lesion")

#### >n calcThresholdMRIaggr ####
# ?calcThresholdMRIaggr
# findMethods("calcThresholdMRIaggr",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, hemisphere = "both", rm.CSF = FALSE, 
#                     threshold = 1:10, decreasing = FALSE, W = NULL, 
#                     GRalgo = FALSE, seed = NULL, name = "Th", trace = TRUE, 
#                     update.object = FALSE, overwrite = FALSE) 

#### test
data(MRIaggr.Pat1_red, package="MRIaggr")

calcW(MRIaggr.Pat1_red,spatial_res=c(1.875,1.875,6),range=6,
      update.object=TRUE,overwrite=TRUE)

calcThresholdMRIaggr(MRIaggr.Pat1_red,param=c("TTP_t0","MTT_t0"),threshold=1:10,
                 rm.CSF=TRUE,hemisphere="lesion",name_newparam=c("TTP.GR_t0","MTT.GR_t0"),
                 GRalgo=TRUE,seed=c("MASK_T2_FLAIR_t2","MASK_DWI_t0"),
                 update.object=TRUE,overwrite=TRUE)


#### example 
## load an MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

#### 1- MRIaggr method 
## raw parameter
multiplot(MRIaggr.Pat1_red,param="TTP_t0",legend=FALSE,main="TTP_t0 - slice ",
             palette=rainbow(10),breaks=seq(0,10)-10^{-10})

## thresholded parameter
calcThresholdMRIaggr(MRIaggr.Pat1_red,param=c("TTP_t0","MTT_t0"),threshold=1:10,
                 name_newparam=c("TTP.th_t0","MTT.th_t0"),
                 update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="TTP.th_t0",main="TTP.th_t0 - slice",
             legend=FALSE,palette=rainbow(10),breaks=(0:10)-10^{-10})

## 1st correction
calcThresholdMRIaggr(MRIaggr.Pat1_red,param=c("TTP_t0","MTT_t0"),threshold=1:10,
                 rm.CSF=TRUE,hemisphere="lesion",name_newparam=c("TTP.red_t0","MTT.red_t0"),
                 update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="TTP.red_t0",main="TTP.red_t0 - slice",
             legend=FALSE,palette=rainbow(10),breaks=(0:10)-10^{-10})

## 2nd correction
calcThresholdMRIaggr(MRIaggr.Pat1_red,param=c("TTP_t0","MTT_t0"),threshold=1:10,
                 rm.CSF=TRUE,hemisphere="lesion",name_newparam=c("TTP.GR_t0","MTT.GR_t0"),
                 GRalgo=TRUE,seed=c("MASK_T2_FLAIR_t2","MASK_DWI_t0"),W.range=sqrt(2),
                 update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="TTP.GR_t0",main="TTP.GR_t0 - slice",
             legend=FALSE,palette=rainbow(10),breaks=(0:10)-10^{-10})

#### >o calcTissueType ####
# ?calcTissueType
# findMethods("calcTissueType",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, niter = 100, nnei = 6, 
#                     beta = "auto", sub = TRUE, trace = TRUE, name_newparam = c("CSF","GM", "WM"),
#                     update.object = FALSE, overwrite = FALSE) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## perform segmentation  (call mritc)
calcTissueType(MRIaggr.Pat1_red,param="T1_t0",update.object=TRUE,overwrite=TRUE,
               niter=100)

## display
multiplot(MRIaggr.Pat1_red,num=1,
             param=c("CSF","WM","GM"),legend=FALSE,
             palette="rgb")

#### >p calcW ####
# ?calcW
# findMethods("calcW",classes="MRIaggr")$MRIaggr
# .local <- function (object, range, spatial_res = c(1, 1, 1), 
#                     num = NULL, hemisphere = "both", subset = NULL, 
#                     upper = TRUE, format = "dgCMatrix", row.norm = FALSE, 
#                     trace = TRUE, update.object = FALSE, overwrite = FALSE) 

#### test

#### example 
## load data and build MRIaggr
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute W (regular lattice)
W <- calcW(MRIaggr.Pat1_red,range=sqrt(2),upper=NULL,num=1:3,hemisphere="lesion")
table(rowSums(W>0))

## compute W (irregular lattice)
W <- calcW(MRIaggr.Pat1_red,range=sqrt(2*1.875^2),upper=NULL,num=1:3,hemisphere="lesion",
           spatial_res=c(1.875,1.875,6))
table(rowSums(W>0))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 4- plot ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### >a boxplotMask ####
# ?boxplotMask
# findMethods("boxplotMask",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, mask, num = NULL, hemisphere = "both", 
#                     norm_mu = FALSE, norm_sigma = FALSE, scale = TRUE, as.logical = FALSE, 
#                     window = FALSE, ylim = NULL, col = c("white", "purple"), 
#                     main = NULL, mgp = c(2, 0.5, 0), x.legend = "topright", 
#                     y.legend = NULL, cex.legend = 0.8, filename = "auto", 
#                     width = 1000, height = 700, path = NULL, unit = "px", 
#                     res = NA) 
  
#### test
boxplotMask(MRIaggr.Pat1_red,param=c("DWI_t0","TTP_t0","MTT_t0"),mask="MASK_T2_FLAIR_t2",
            as.logical=TRUE,x.legend="bottomleft")

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## display
boxplotMask(MRIaggr.Pat1_red,param=c("DWI_t0","TTP_t0","MTT_t0"),mask="MASK_T2_FLAIR_t2",
            as.logical=TRUE)

#### >b heatmapMRIaggr ####
# ?heatmapMRIaggr
# findMethods("heatmapMRIaggr",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, num = NULL, hemisphere = "both", 
#                     scale = TRUE, method = "pearson", points.values = TRUE, 
#                     type = "image", digit = 3, breaks = NULL, window = FALSE, 
#                     col = cm.colors(256), main = NULL, mgp = c(2, 0.5, 0), 
#                     mar = c(4, 4, 1, 6), las = 1, cex.axis = 1, filename = "auto", 
#                     width = 1000, height = 700, path = NULL, unit = "px", 
#                     res = NA) 

#### test
# kendall : long !!!
heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image",cex=0.75,method="spearman",
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red",package="MRIaggr")

## pearson
heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image",cex=0.75,
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))    
## spearman
heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image",cex=0.75,method="spearman",
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))  


## spearman with legend
heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image.plot",cex=0.75,method="spearman",
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))  


#### >c multiplot ####
# ?multiplot
# findMethods("multiplot",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, num = NULL, index1 = NULL, 
#                     index2 = NULL, index3 = NULL, midplane = FALSE, 
#                     slice_var = "k", hemisphere = "both", norm_mu = FALSE, 
#                     norm_sigma = FALSE, as.logical = FALSE, breaks = 50, 
#                     type.breaks = "range", palette = "terrain.colors", col = NULL, 
#                     pch = 15, cex = 1, col.NA = "lightyellow", pch.NA = 8, 
#                     col.midplane = "red", xlim = NULL, ylim = NULL, axes = TRUE, 
#                     window = FALSE, legend = TRUE, mfrow = NULL, mar = rep(1.5,4), mgp = c(2, 0.5, 0), pty = NULL, asp = 1, bg = "lightblue", 
#                     xlab = "", ylab = "", main = paste(" slice ", slice_var,sep = ""), num.main = TRUE, cex.main = 1.5, 
#                     quantiles.legend = TRUE, digit.legend = 3, cex.legend = 1.5, 
#                     mar.legend = c(1.5, 7, 1.5, 1.5), main.legend = param, 
#                     filename = "multiplot", width = 1000, height = 700, 
#                     path = NULL, unit = "px", res = NA) 

#### test

# gestion de la fenetre graphique
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
legend=F,window="png")
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
legend=T,window="png")
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
legend=NULL,window="png")


#### example 
## load an MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

# display 3 slices 
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3)

# display 3 slices with no axes and white background
multiplot(MRIaggr.Pat1_red,param="DWI_t0",
             num=1:3,axes=FALSE,bg="white")

# remove the legend
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3,legend=FALSE)

## display an set of points
# using a binary parameter stored in the object
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3,index1=list(coords="MASK_DWI_t0")
)

# customize the display of the points
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3,index1=list(coords="MASK_DWI_t0",col="pink",pch=14)
)

# display only the edges of the set
multiplot(MRIaggr.Pat1_red,param="DWI_t0",num=3, legend=FALSE,  
             index1=list(coords="MASK_DWI_t0",outline=TRUE)
)


# specify the index of points using coordinates
coordsIndex <- data.frame(i=c(40,60),j=c(80,100),k=c(3,3))

multiplot(MRIaggr.Pat1_red,param="DWI_t0",num=3,legend=FALSE,       
             index2=list(coords=coordsIndex,col="black",pch=15,cex=4)
)


# various possibilities for the display
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
             legend=FALSE,window=FALSE)
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
             legend=TRUE,window=FALSE)
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
             legend=NULL,window=FALSE)


#### >d outlineMRIaggr ####
# ?outlineMRIaggr
# findMethods("outlineMRIaggr",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, num = NULL, hemisphere = "both", 
#                     xlim = NULL, ylim = NULL, palette.multiplot = "terrain.colors", 
#                     col.multiplot = NULL, breaks.multiplot = 25, fill = TRUE, 
#                     n = 50, sequential = FALSE, min_dist = 1, col = c("blue", 
#                                                                       "red", "grey"), pch = 20, cex = c(0.75, 1, 0.75), 
#                     trace = TRUE, update.object = FALSE, overwrite = FALSE, 
#                     name_newparam = "userMask") 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## outline the area of interest
res <- outlineMRIaggr(MRIaggr.Pat1_red,param="DWI_t0",
                  num=3:5,sequential=TRUE,overwrite=TRUE,update.object=TRUE)

multiplot(MRIaggr.Pat1_red,param="userMask",              
             num=3:5)


## outline an edge of interest 
res <- outlineMRIaggr(MRIaggr.Pat1_red,param="DWI_t0",index1="MASK_DWI_t0",
                  fill=FALSE,num=3:5,sequential=TRUE,overwrite=TRUE,update.object=TRUE)

multiplot(MRIaggr.Pat1_red,param="userMask",              
             num=3:5)

## define an new area as the union of the outlined area and the initial lesion mask
res <- outlineMRIaggr(MRIaggr.Pat1_red,param="DWI_t0",index1="MASK_DWI_t0",
                  operator_index1="union",num=3,sequential=TRUE,overwrite=TRUE,update.object=TRUE)

multiplot(MRIaggr.Pat1_red,param="userMask",              
             num=3)

## define an new area as the intersection of the outlined area and the initial lesion mask
res <- outlineMRIaggr(MRIaggr.Pat1_red,param="DWI_t0",index1="MASK_DWI_t0",
                  operator_index1="intersection",num=3,sequential=TRUE,overwrite=TRUE,update.object=TRUE)

multiplot(MRIaggr.Pat1_red,param="userMask",              
             num=3)

## define an new area as the difference between the outlined area and the initial lesion mask
res <- outlineMRIaggr(MRIaggr.Pat1_red,param="DWI_t0",index1="MASK_DWI_t0",
                  operator_index1="difference",num=3,sequential=TRUE,overwrite=TRUE,update.object=TRUE)

multiplot(MRIaggr.Pat1_red,param="userMask",              
             num=3)

#### >e plotDistClass ####
# ?plotDistClass
# findMethods("plotDistClass",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, class, num = NULL, hemisphere = "both", 
#                     norm_mu = FALSE, norm_sigma = FALSE, bw.adjust = 1, kernel = "gaussian", 
#                     from = NULL, to = NULL, ylim = NULL, window = FALSE, 
#                     col = 1:6, main = NULL, mgp = c(2, 0.5, 0), type = "l", 
#                     pch = 20, lwd = 1, x.legend = "topright", y.legend = NULL, 
#                     cex.legend = 0.8, filename = "auto", width = 1000, height = 700, 
#                     path = NULL, unit = "px", res = NA) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## display
plotDistClass(MRIaggr.Pat1_red,param="DWI_t0",
              class=c("MASK_T2_FLAIR_t2"))

# specify the smoothing bandwidth
plotDistClass(MRIaggr.Pat1_red,param="DWI_t0",bw.adjust=1,
              class=c("MASK_T2_FLAIR_t2"))

# specify the scale
plotDistClass(MRIaggr.Pat1_red,param="DWI_t0",bw.adjust=1,
              from=200,to=300,class=c("MASK_T2_FLAIR_t2"))

# usee several classes
plotDistClass(MRIaggr.Pat1_red,param="TTP_t0",bw.adjust=2,
              class=c("CSF","WM","GM","MASK_T2_FLAIR_t2"))

plotDistClass(MRIaggr.Pat1_red,param="DWI_t0",bw.adjust=2,
              class=c("CSF","WM","GM","MASK_T2_FLAIR_t2"))

#### >f plotLesion3D ####
# ?plotLesion3D
# findMethods("plotLesion3D",classes="MRIaggr")$MRIaggr
# .local <- function (object, mask, edge = FALSE, Neighborhood = "3D_N6", 
#                     as.logical = FALSE, spatial_res = c(1, 1, 1), xlim = NULL, 
#                     ylim = NULL, zlim = NULL, type.plot = "shapelist3d", 
#                     px_max = 10000, radius = 1, type = "s", col = "red", 
#                     col.edge = "black") 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

if(require(rgl)){
  # global view
  plotLesion3D(MRIaggr.Pat1_red,mask="MASK_T2_FLAIR_t2",spatial_res=c(1.875,1.875,6),
               as.logical=TRUE)
  
  # by slice
  plotLesion3D(MRIaggr.Pat1_red,mask="MASK_T2_FLAIR_t2",spatial_res=c(1.875,1.875,6),
               type.plot="plot3d",
               as.logical=TRUE)
}

#### >g plotTableLesion ####
# ?plotTableLesion
# findMethods("plotTableLesion",classes="MRIaggr")$MRIaggr
# .local <- function (object, mask, num = NULL, type = "matplot", 
#                     window = FALSE, col = 1:5, lty = 1:5, lwd = 1, main = NULL, 
#                     mgp = c(2, 0.5, 0), mar = rep(3, 4), cex.legend = 1, 
#                     cex.main = 1, cex.axis = 1, cex.lab = 1, filename = "auto", 
#                     width = 1000, height = 700, path = NULL, unit = "px", 
#                     res = NA) 

#### test

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## matplot display of the lesion
plotTableLesion(MRIaggr.Pat1_red,num=1:3,type="matplot",
                mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"))

## evolution display of the lesion
plotTableLesion(MRIaggr.Pat1_red,num=1:3,type="evolution",
                mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"))

#### >h pointsHemisphere ####
# ?pointsHemisphere
# findMethods("pointsHemisphere",classes="MRIaggr")$MRIaggr
# .local <- function (object, col = "red", lwd = 2, lty = 1) 

#### test
# pointsHemisphere(MRIaggr.Pat1)

#### example 


#### >i summary ####
# ?"summary,MRIaggr-method"
# findMethods("summary",classes="MRIaggr")$MRIaggr
# .local <- function (object, param = FALSE, clinic = FALSE, 
#                     descStats = FALSE, history = FALSE) 

#### test
summary(MRIaggr.Pat1_red,param=TRUE)

summary(MRIaggr.Pat1_red,param=TRUE,clinic=TRUE)

summary(MRIaggr.Pat1_red,param=TRUE,clinic=TRUE,
        descStats = TRUE, history = TRUE)

#### example 
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

summary(MRIaggr.Pat1_red)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 5- const ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#### >a constCompressMRIaggr ####
# ?constCompressMRIaggr
# findMethods("constCompressMRIaggr",classes="MRIaggr")$MRIaggr
# .local <- function (object, factor, param=NULL, mask = NULL, threshold = 0.49, 
#                     trace = FALSE) 

#### test

#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compress the MRIaggr object
MRIaggr.compressed <- constCompressMRIaggr(MRIaggr.Pat1_red,factor=2,
                                           param=c("DWI_t0","T2_FLAIR_t2","MASK_T2_FLAIR_t2"),
                                           mask="MASK_T2_FLAIR_t2") 

## display
par(mfrow=c(2,4),mar=rep(1.75,4),mgp=c(2,0.75,0))
multiplot(MRIaggr.Pat1_red,param="DWI_t0",window=NULL,breaks=seq(0,350,1),
             midplane=TRUE,main="before - slice ")
multiplot(MRIaggr.compressed,param="DWI_t0",window=NULL,breaks=seq(0,350,1),
             midplane=TRUE,main="after - slice ")

multiplot(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2",main="before - slice ")
multiplot(MRIaggr.compressed,param="MASK_T2_FLAIR_t2",main="after - slice ")

#### >b constReduceMRIaggr ####
# ?constReduceMRIaggr
# findMethods("constReduceMRIaggr",classes="MRIaggr")$MRIaggr
# .local <- function (object, mask, as.logical = FALSE, keep.index = TRUE) 

#### test

#### example
## load NIFTI files and convert them to MRIaggr
path <- system.file(file.path("nifti"),package = "MRIaggr")
ls.array <- list(readMRI(file=file.path(path,"T2_GRE_t0"),format="nifti"))
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param="T2_GRE_t0")

## create the cerebral mask
res <- calcBrainMask(MRIaggr.Pat1,param="T2_GRE_t0",type="kmeans",
                     kmeans.n_groups=2:4,
                     update.object=TRUE,overwrite=TRUE)

res <- calcSmoothMask(MRIaggr.Pat1,update.object=TRUE,overwrite=TRUE)

## display
multiplot(MRIaggr.Pat1,param="mask",legend=FALSE)  				  

## construct the reduced object
MRIaggr.Pat1_red <- constReduceMRIaggr(MRIaggr.Pat1,mask="mask")

## display
par(mfrow=c(2,4),mar=rep(1.75,4),mgp=c(2,0.75,0))
multiplot(MRIaggr.Pat1,num=c(1,3,6),param="T2_GRE_t0",
             window=NULL,breaks=seq(0,300,1),legend=FALSE)
multiplot(MRIaggr.Pat1_red,num=c(1,3,6),param="T2_GRE_t0",
             window=NULL,breaks=seq(0,300,1),legend=FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 6- init ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#### 6.a initParameter - private ####
# ?initParameter
# findMethods("initParameter",classes="MRIaggr")$MRIaggr
# .local <- function (object, param, test = TRUE, init = FALSE, 
#                     accept.coords = TRUE, accept.mask = TRUE, accept.index = TRUE, 
#                     arg_name = "param", long_name = "parameters", method = "initParameter")
  

#### test
data("MRIaggr.Pat1_red",package="MRIaggr")
res <- MRIaggr:::initParameter(MRIaggr.Pat1_red,param="DWI_t0")

#### example

#### 6.b initNum - private #### 
# ?initNum
# findMethods("initNum",classes="MRIaggr")$MRIaggr
# .local <- function (object, num, test = TRUE, init = TRUE, 
#                     slice_var = "k", method = "initNum") 

#### test
MRIaggr:::initNum(MRIaggr.Pat1_red,num=NULL,init=FALSE)
MRIaggr:::initNum(MRIaggr.Pat1_red,num=NULL)
  
MRIaggr:::initNum(MRIaggr.Pat1_red,num=1:3,init=FALSE)
MRIaggr:::initNum(MRIaggr.Pat1_red,num=1:3,init=FALSE,slice_var="i")
MRIaggr:::initNum(MRIaggr.Pat1_red,num=-1:3,init=FALSE,test=FALSE)

#### example




