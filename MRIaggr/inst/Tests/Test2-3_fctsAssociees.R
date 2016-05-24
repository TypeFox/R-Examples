#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Fichier test 2.3 : Test des fonctions associees 
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require(MRIaggr)

options(error=function() traceback(2)) 
options(max.print=10000)

#### 0- chargement ####

#### 1- Const ####
#### >a constCarto3D - private  ####
# ?constCarto3D
# MRIaggr:::constCarto3D
# function (array, identifier, parameter, default_value = NULL, 
#           pos_default_value = c(1, 1, 1), rm.array = FALSE) 

#### test
path <- system.file(file.path("nifti"),package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")

Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,
                                    identifier="Pat1",param="TTP_t0",default_value=NULL,
                                    pos_default_value=cbind(1:10,1,1))

Carto3D.Pat1_TTP_t0@default_value
Carto3D.Pat1_TTP_t0@identifier
Carto3D.Pat1_TTP_t0@parameter
nifti.Pat1_TTP_t0[,1,1,1]

Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,
                                    identifier="Pat1",param="TTP_t0",default_value=NA)

Carto3D.Pat1_TTP_t0@default_value

Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,
                                    identifier="Pat1",param="TTP_t0",default_value=NULL,
                                    pos_default_value=cbind(1:10,1,1), rm.array = TRUE)

"nifti.Pat1_TTP_t0" %in% ls()

# example 
## load NIFTI files
path <- system.file(file.path("nifti"),package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")

## convert them to Carto3D
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")
class(Carto3D.Pat1_TTP_t0)

#### >b constMRIaggr ####
# ?constMRIaggr
# constMRIaggr
# function (ls.array, identifier, parameters, default_value = NULL, 
#           pos_default_value = c(1, 1, 1), tol = 10^{-10}, trace = TRUE, rm.ls.array = FALSE) 
  

#### test
path <- system.file(file.path("nifti"),package = "MRIaggr")

param <- c("DWI_t0.nii","MASK_DWI_t0.nii","MTT_t0.nii","TTP_t0.nii","T1_t0.nii","T2_GRE_t0.nii",
           "MTT_t1.nii","TTP_t1.nii","T2_FLAIR_t2.nii","MASK_T2_FLAIR_t2.nii")

ls.array <- list()
for(iter_param in 1:length(param)){
  ls.array[[iter_param]] <- readMRI(file=file.path(path,param[iter_param]),format="nifti")
}

MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=param,default_value=5:15)
MRIaggr.Pat1@default_value
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=param,default_value=5:15,trace=FALSE)
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=param,pos_default_value=cbind(5:15,1,1),trace=FALSE)
MRIaggr.Pat1@default_value
table(ls.array[[1]][cbind(5:15,1,1,1)])
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=param,pos_default_value=cbind(5:15,1,1),rm.ls.array=TRUE)

#### example 
#### 1- 1st method
## load NIFT files
path <- system.file(file.path("nifti"),package = "MRIaggr")

nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")
nifti.Pat1_DWI_t0 <- readMRI(file=file.path(path,"DWI_t0"),format="nifti")
nifti.Pat1_MASK_DWI_t0 <- readMRI(file=file.path(path,"MASK_DWI_t0"),format="nifti")
nifti.Pat1_MASK_T2_FLAIR_t2 <- readMRI(file=file.path(path,"MASK_T2_FLAIR_t2"),format="nifti")

## convert them to MRIaggr
MRIaggr.Pat1 <- constMRIaggr(list(nifti.Pat1_TTP_t0,nifti.Pat1_DWI_t0,
                                  nifti.Pat1_MASK_DWI_t0,nifti.Pat1_MASK_T2_FLAIR_t2),
                             identifier="Pat1",param=c("TTP_t0","DWI_t0","MASK_DWI_t0","MASK_T2_FLAIR_t2"))


#### 2- 2nd method
## load NIFTI files
param <- c("DWI_t0.nii","MASK_DWI_t0.nii","MTT_t0.nii","TTP_t0.nii","T1_t0.nii","T2_GRE_t0.nii",
           "MTT_t1.nii","TTP_t1.nii","T2_FLAIR_t2.nii","MASK_T2_FLAIR_t2.nii")

ls.array <- list()
for(iter_param in 1:length(param)){
  ls.array[[iter_param]] <- readMRI(file=file.path(path,param[iter_param]),format="nifti")
}

## convert them to MRIaggr
param <- gsub(".nii","",param)

MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=param)

#### >c constSweave ####
# ?constSweave
# constSweave
# function (dir, identifier = NULL, param = NULL, table = NULL, 
#           extra_text = NULL, subsection = NULL, index_subsection = NULL, 
#           subsubsection = NULL, index_subsubsection = NULL, legend = NULL, 
#           trace = FALSE, width = list(0.9, 0.9, 0.9), trim = list(c(0,0, 0, 0), c(0, 0, 0, 0), c(0, 160, 0, 0)), width.legend = 0.35, 
#           trim.legend = c(0, 0, 0, 0), title = "", date = "", author = "") 

#### test : cas special - traite a part

#### example : cas special - traite a part

#### 2- Conversion ####

#### >a array2df ####
# ?array2df
# array2df
# function (array, coords = NULL, param = "res", names_coords = letters[9:(8 +ncol(coords))], na.rm = TRUE) 
 
#### test
set.seed(10)
n <- 4
Y <- rnorm(n^2)

Yarray <- array(Y,dim=c(4,4))

res1 <- array2df(array=Yarray,coords=expand.grid(1:n+0.5,1:n+0.5))
res1-cbind(expand.grid(1:n+0.5,1:n+0.5),Y)
res2 <- array2df(array=Yarray)
res2-cbind(expand.grid(1:n,1:n),Y)

Yarray <- array(NA,dim=c(10,10))
grid <- expand.grid(2*(1:n),2*(1:n))
Yarray[ apply(which(is.na(Yarray),arr.ind=TRUE),1,function(x){(x[1] %in% grid[,1])*(x[2] %in% grid[,2])})  == 1] <- Y
res4 <- array2df(array=array(c(Yarray,Yarray,Yarray),dim=c(10,10,3)))
res4[res4$k==1,]-res4[res4$k==2,]

#### example 
## load a NIFTI file (array format)
path <- system.file(file.path("nifti"),package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")
dim(nifti.Pat1_TTP_t0)

## conversion to data frame format
res128 <- array2df(array=nifti.Pat1_TTP_t0,name_newparam="TTP_t0")
dim(res128)
head(res128)

## conversion to data frame format with specific coordinates
res256 <- array2df(array=nifti.Pat1_TTP_t0,name_newparam="TTP_t0",
                   coords=expand.grid(129:206,129:228,1:3,1))
dim(res256)
head(res256)


#### >b Carto3D2MRI - private ####
# ?Carto3D2MRIaggr
# MRIaggr:::Carto3D2MRIaggr
# function (ls.Carto3D, rm.Carto3D = FALSE, tol = 10^{-10}, num = NULL, trace = TRUE) 
  
#### test
path <- system.file(file.path("nifti"),package = "MRIaggr")

Pat1.TTP.t0.nifti <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")
Pat1.DWI.t0.nifti <- readMRI(file=file.path(path,"DWI_t0"),format="nifti")
Pat1.MASK_DWI.t0.nifti <- readMRI(file=file.path(path,"MASK_DWI_t0"),format="nifti")
Pat1.MASK_T2_FLAIR.t2.nifti <- readMRI(file=file.path(path,"MASK_T2_FLAIR_t2"),format="nifti")

Pat1.TTP.t0.Carto3D <- constCarto3D(Pat1.TTP.t0.nifti,
                                    identifier="Pat1",param="TTP_t0",default_value=NA)
Pat1.DWI.t0.Carto3D <- constCarto3D(Pat1.DWI.t0.nifti,
                                    identifier="Pat1",param="DWI_t0",default_value=NA)
Pat1.MASK_DWI.t0.Carto3D <- constCarto3D(Pat1.MASK_DWI.t0.nifti,
                                         identifier="Pat1",param="MASK_DWI_t0",default_value=NA)
Pat1.MASK_T2_FLAIR.t2.Carto3D <- constCarto3D(Pat1.MASK_T2_FLAIR.t2.nifti,
                                              identifier="Pat1",param="MASK_T2_t2",default_value=NA)

MRIaggr.Pat1 <- Carto3D2MRIaggr(list(Pat1.TTP.t0.Carto3D,
                                     Pat1.DWI.t0.Carto3D,
                                     Pat1.MASK_DWI.t0.Carto3D,
                                     Pat1.MASK_T2_FLAIR.t2.Carto3D),
                                num=1:3
)

MRIaggr.Pat1 <- Carto3D2MRIaggr(list(Pat1.TTP.t0.Carto3D,
                                     Pat1.DWI.t0.Carto3D,
                                     Pat1.MASK_DWI.t0.Carto3D,
                                     Pat1.MASK_T2_FLAIR.t2.Carto3D),
                                num=1:3,rm.Carto3D=TRUE
)

# example
## load NIFTI files 
path <- system.file(file.path("nifti"),package = "MRIaggr")

Pat1.TTP.t0.nifti <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")
Pat1.DWI.t0.nifti <- readMRI(file=file.path(path,"DWI_t0"),format="nifti")
Pat1.MASK_DWI.t0.nifti <- readMRI(file=file.path(path,"MASK_DWI_t0"),format="nifti")
Pat1.MASK_T2_FLAIR.t2.nifti <- readMRI(file=file.path(path,"MASK_T2_FLAIR_t2"),format="nifti")

## convert them to Carto3D
Pat1.TTP.t0.Carto3D <- constCarto3D(Pat1.TTP.t0.nifti,
                                    identifier="Pat1",param="TTP_t0",default_value=NA)
Pat1.DWI.t0.Carto3D <- constCarto3D(Pat1.DWI.t0.nifti,
                                    identifier="Pat1",param="DWI_t0",default_value=NA)
Pat1.MASK_DWI.t0.Carto3D <- constCarto3D(Pat1.MASK_DWI.t0.nifti,
                                         identifier="Pat1",param="MASK_DWI_t0",default_value=NA)
Pat1.MASK_T2_FLAIR.t2.Carto3D <- constCarto3D(Pat1.MASK_T2_FLAIR.t2.nifti,
                                              identifier="Pat1",param="MASK_T2_t2",default_value=NA)

## convert Carto3D to MRIaggr  							 
MRIaggr.Pat1 <- Carto3D2MRIaggr(list(Pat1.TTP.t0.Carto3D,
                                     Pat1.DWI.t0.Carto3D,
                                     Pat1.MASK_DWI.t0.Carto3D,
                                     Pat1.MASK_T2_FLAIR.t2.Carto3D)
)


#### >c df2array ####
# ?df2array
# df2array
# function (data, coords, format = "any", default_value = NA, range.coords = NULL) 
  
#### test
set.seed(10)
n <- 4
Y <- rnorm(n^2)
res1 <- df2array(contrast=Y,coords=expand.grid(1:n+0.5,1:n+0.5))
res2 <- df2array(contrast=Y,coords=expand.grid(1:n,1:n),format="matrix")
res3 <- df2array(contrast=Y,coords=expand.grid(2*(1:n),2*(1:n)))
res4 <- df2array(contrast=cbind(Y,Y,Y),coords=expand.grid(2*(1:n),2*(1:n)),
                 range.coords=c(10,10))

par(mfrow=c(2,2),mar=rep(2,4),mgp=c(1.5,0.5,0))
fields::image.plot(unique(res1$coords[,1]),unique(res1$coords[,2]),res1$contrast[[1]],xlab="",ylab="")
fields::image.plot(unique(res2$coords[,1]),unique(res2$coords[,2]),res2$contrast,xlab="",ylab="")
fields::image.plot(res3$contrast[[1]])
fields::image.plot(res4$contrast[[2]])


data("MRIaggr.Pat1_red",package="MRIaggr")
carto <- selectContrast(MRIaggr.Pat1_red,param="DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)
array.DWI_t0 <- df2array(carto,coords=coords,default_value=1000)$contrast[[1]]
fields::image.plot(array.DWI_t0[,,1])
array.DWI_t0 <- df2array(carto,coords=coords,default_value=1000, 
                         range.coords = c(256,256,6))$contrast[[1]]
fields::image.plot(array.DWI_t0[,,1])


#### example
#### 1- with simulated data
## simulate
set.seed(10)
n <- 4
Y <- rnorm(n^2)

## conversion
res1 <- df2array(contrast=Y,coords=expand.grid(1:n+0.5,1:n+0.5))
res2 <- df2array(contrast=Y,coords=expand.grid(1:n,1:n),format="matrix")
res3 <- df2array(contrast=Y,coords=expand.grid(2*(1:n),2*(1:n)))
res4 <- df2array(contrast=cbind(Y,Y,Y),coords=expand.grid(2*(1:n),2*(1:n)),
                 range.coords=c(10,10))

## display
par(mfrow=c(2,2),mar=rep(2,4),mgp=c(1.5,0.5,0))
fields::image.plot(unique(res1$coords[,1]),unique(res1$coords[,2]),res1$contrast[[1]],xlab="",ylab="")
fields::image.plot(unique(res2$coords[,1]),unique(res2$coords[,2]),res2$contrast,xlab="",ylab="")
fields::image.plot(res3$contrast[[1]])
fields::image.plot(res4$contrast[[2]])

#### 2- with MRIaggr data
## load an MRIaggr object
data("MRIaggr.Pat1_red",package="MRIaggr")
carto <- selectContrast(MRIaggr.Pat1_red,param="DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)

## converion 1
array.DWI_t0 <- df2array(carto,coords=coords,default_value=1000)$contrast[[1]]

# display
fields::image.plot(array.DWI_t0[,,1])

## conversion 2
array.DWI_t0 <- df2array(carto,coords=coords,default_value=1000, 
                         range.coords = c(256,256,6))$contrast[[1]]

# display
fields::image.plot(array.DWI_t0[,,1])

  
#### 3- calc ####
#### >a calcAUPRC ####
# ?calcAUPRC
# calcAUPRC
# function (x, y, subdivisions = 10000, performance = NULL)

#### test
n0 <- 1000
n1 <-1000
X <- c(rnorm(n0,0),rnorm(n1,2))
Y <- c(rep(0,n0),rep(1,n1))

calcAUPRC(X,Y,subdivisions=10000)
calcAUPRC(X,Y,subdivisions=1000)
calcAUPRC(X,Y,subdivisions=100)
# calcAUPRC(X,Y,subdivisions=50)
# calcAUPRC(X,Y,subdivisions=25)
# calcAUPRC(X,Y,subdivisions=20)

perfXY <- ROCR::performance(ROCR::prediction(X,Y),x.measure="rec",measure="prec")
# perfXY <- performance(prediction(X,Y),x.measure="prec",measure="err")
calcAUPRC(performance=perfXY,subdivisions=10000)


#### example
#### 1- with MRIaggr data 
## load an MRIaggr object
data(MRIaggr.Pat1_red,package="MRIaggr")

## select parameter and binary outcome
cartoT2 <- selectContrast(MRIaggr.Pat1_red,param="T2_FLAIR_t2")
cartoMASK <- selectContrast(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2")

## compute AUPRC
T2.AUPRC <- calcAUPRC(x=cartoT2,y=cartoMASK)

## compute AUC
# if(require(pROC)){
# T2.AUC <- auc(roc(cartoMASK ~ cartoT2))
# }

## display
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=1,
             index1=list(coords="MASK_T2_FLAIR_t2",outline=TRUE)
)

#### 2- with simulated data 
n0 <- 1000
n1 <- c(10,100,1000)
for(iter_n in 1:length(n1)){
  X <- c(rnorm(n0,0),rnorm(n1[iter_n],2))
  Y <- c(rep(0,n0),rep(1,n1[iter_n]))
  print(calcAUPRC(X,Y))
}

## alternative way using a performance object
perfXY <- ROCR::performance(ROCR::prediction(X,Y),x.measure="rec",measure="prec")
calcAUPRC(performance=perfXY,subdivisions=10000)

#### >b calcFilter ####
# ?calcFilter
# findMethods(MRIaggr:::calcFilter,classes="array")$array
# .local <- function (object, filter, w_contrast = FALSE, na.rm = FALSE) 

#### test

n <- 100

#### test 1
M <- matrix(rnorm(n),nrow=sqrt(n),ncol=sqrt(n))
M[1:5,1:5] <- rnorm(25,mean=10,sd=2)
M[8,8] <- -7
A <- array(c(rep(1,length(M)),M,rep(10,length(M))),dim=c(10,10,3))

#### test 2
M <- matrix(rnorm(n),nrow=sqrt(n),ncol=sqrt(n))
M[] <- 1
M[8,8] <- -7
A <- array(c(1,M,1),dim=c(10,10,3))
A[5,5,1:3] <- rnorm(3,-5)
A[1:3,1,2] <- rnorm(3,5)
A[4,7:9,2] <- rnorm(3,5)

filtreG3_F.C <- calcFilter(M,filter="2D_G3",na.rm=F)$res
filtreG3_T.C <- calcFilter(M,filter="2D_G3",na.rm=T)$res

filtreG3med_F.C <- calcFilter(M,filter="2D_M3",na.rm=F)$res
filtreG3w_F.C <- calcFilter(M,filter="2D_G3",na.rm=F,w_contrast=T)$res

filtre3DG3med_F.C <- calcFilter(A,filter="3D_M3",na.rm=F)$res
filtre3DG3_F.C <- calcFilter(A,filter="3D_G3",na.rm=F)$res
filtre3DG3w_F.C <- calcFilter(A,filter="3D_G3",na.rm=F,w_contrast=T)$res

if(try(require(fields))){
  breaks <- seq(min(M),max(M),length.out=200)
  n.breaks <- length(breaks)
  
  par(mfrow=c(2,3),mar=rep(2,4),mgp=c(1.5,0.5,0))
  fields::image.plot(M,breaks=breaks, col=terrain.colors(n.breaks-1),main="initial data")
  fields::image.plot(filtreG3_F.C,breaks=breaks, col=terrain.colors(n.breaks-1),main="gaussian filter (bdw 3px/na.rm=F/.cpp)",zlim=range(M))
  fields::image.plot(filtreG3_T.C,breaks=breaks, col=terrain.colors(n.breaks-1),main="gaussian filter bandwidth 3px (bdw 3px/na.rm=T/.cpp)",zlim=range(M))
  
  fields::image.plot(filtreG3w_F.C,breaks=breaks, col=terrain.colors(n.breaks-1),main="Wgaussian filter (bdw 3px/na.rm=F/.cpp)",zlim=range(M))
  fields::image.plot(filtreG3med_F.C,breaks=breaks, col=terrain.colors(n.breaks-1),main="median filter (bdw 3px/na.rm=F/.cpp)",zlim=range(M))
  
  
  breaks <- seq(min(A),max(A),length.out=200)
  n.breaks <- length(breaks)
  
  par(mfrow=c(2,2),mar=rep(2,4),mgp=c(1.5,0.5,0))
  fields::image.plot(A[,,2],breaks=breaks, col=terrain.colors(n.breaks-1),main="initial data",zlim=range(A))
  fields::image.plot(filtre3DG3_F.C[,,2],breaks=breaks, col=terrain.colors(n.breaks-1),main="gaussian filter (bdw 3px/na.rm=F/.cpp)",zlim=range(A))
  fields::image.plot(filtre3DG3w_F.C[,,2],breaks=breaks, col=terrain.colors(n.breaks-1),main="Wgaussian filter (bdw 3px/na.rm=F/.cpp)",zlim=range(A))
  fields::image.plot(filtre3DG3med_F.C[,,2],breaks=breaks, col=terrain.colors(n.breaks-1),main="median filter (bdw 3px/na.rm=F/.cpp)",zlim=range(A))
}

#### test 3
M <- matrix(1,nrow=5,ncol=5)
filtre2DN.C <- calcFilter(M,filter="2D_N4")$res
filtre2Dn.C <- calcFilter(M,filter="2D_N4")$res
filtre2DN_narm.C <- calcFilter(M,filter="2D_N4",na.rm=TRUE)$res
filtre2Dn_narm.C <- calcFilter(M,filter="2D_N4",na.rm=TRUE)$res

#### test 4
# data("MRIaggr.Pat1_red", package="MRIaggr") 


#### example
## load a NIFTI file
path <- system.file(file.path("nifti"),package = "MRIaggr")
nifti.Pat1_DWI_t0 <- readMRI(file=file.path(path,"DWI_t0"),format="nifti")

## before filtering
graphics::image(nifti.Pat1_DWI_t0[,,1,1])

## after median filtering
niftiF.Pat1_DWI_t0 <- calcFilter(nifti.Pat1_DWI_t0[,,,1],filter="2D_M3")$res
graphics::image(niftiF.Pat1_DWI_t0[,,1])

#### >c calcGroupsCoords ####
# ?calcGroupsCoords
# calcGroupsCoords
# function (coords, array = NULL, Neighborhood, max_groups = 10000, 
#           trace = TRUE) 

#### test
n <- 10
coords <- rbind(which(matrix(0,nrow=n,ncol=n)==0,arr.ind = TRUE),
                which(matrix(0,nrow=n,ncol=n)==0,arr.ind = TRUE)+100
)
res.Groups <- calcGroupsCoords(coords=coords,
                               Neighborhood="2D_N4")

W <- calcW(object=as.data.frame(coords),range=sqrt(2),row.norm=TRUE)
if(try(require(fields))){
  fields::image.plot(as.matrix(W))
}

Acoords <- df2array(rep(1,200),coords)
res.Groups <- calcGroupsCoords(array=Acoords$contrast[[1]],
                               Neighborhood="2D_N4")

names(res.Groups)
res.Groups$group_size


#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select data
MASK_DWI_t0 <- selectContrast(MRIaggr.Pat1_red,param="MASK_DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)

#### 1- compute spatial groups using coordinates
res3DN18 <- calcGroupsCoords(coords=coords[MASK_DWI_t0==1,],Neighborhood="3D_N18")
res3DN18$group_size

## display the lesion within the brain
multiplot(coords,contrast=MASK_DWI_t0,
             legend=FALSE,xlim=c(30,100),ylim=c(20,110),
             index1=list(coords=coords,outline=TRUE,cex=0.5),
             num.main=FALSE)

## display the lesion spatial groups
multiplot(res3DN18$df.group[,c("i","j","k")],contrast=res3DN18$df.group[,"group"],
             legend=FALSE,xlim=c(30,100),ylim=c(20,100),cex=0.5,
             index1=list(coords=coords,outline=TRUE))

## zoom
multiplot(res3DN18$df.group[,c("i","j","k")],contrast=res3DN18$df.group[,"group"],
             legend=FALSE,xlim=c(30,100),ylim=c(20,100),num=1,cex=0.1,
             index1=list(coords=coords,outline=TRUE))

#### 2-compute spatial groups using an array
A.MASK_DWI_t0 <- df2array(contrast=MASK_DWI_t0,coords)$contrast[[1]]
A.MASK_DWI_t0[A.MASK_DWI_t0==FALSE] <- NA

## display
graphics::image(A.MASK_DWI_t0[,,3])

## computation of the spatial groups
res3DN18.bis <- calcGroupsCoords(array=A.MASK_DWI_t0,Neighborhood="3D_N18")

res3DN18$group_size-res3DN18.bis$group_size # same result

#### >d calcGroupsW ####
# ?calcGroupsW
# calcGroupsW
# function (W, subset = NULL, max_groups = 10000) 

#### test
if(require(spam) && require(Matrix)){
  n <- 10
  coords <- rbind(which(matrix(0,nrow=n,ncol=n)==0,arr.ind = TRUE),
                  which(matrix(0,nrow=n,ncol=n)==0,arr.ind = TRUE)+100
  )
  
  W <- calcW(object=as.data.frame(coords),range=sqrt(2),row.norm=TRUE)
  if(try(require(fields))){
    fields::image.plot(as.matrix(W))
  }
  
  res.Groups <- calcGroupsW(W)
  
  names(res.Groups)
  res.Groups$group_size
}


#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select data
MASK_DWI_t0 <- selectContrast(MRIaggr.Pat1_red,param="MASK_DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)

## select compute W
W <- calcW(object=as.data.frame(coords[MASK_DWI_t0==1,]),range=sqrt(2),row.norm=TRUE)

## find spatial groups
res.Groups <- calcGroupsW(W)
res.Groups$group_size

## display
multiplot(coords[MASK_DWI_t0==1,],contrast=res.Groups$group,
             legend=FALSE,xlim=c(30,100),ylim=c(20,100),cex=0.5,
             palette=rainbow(10)[-1],index1=list(coords=coords,outline=TRUE))


#### >e calcThreshold #### 
# ?calcThreshold
# calcThreshold
# function (data, param, hemisphere = NULL, rm.CSF = FALSE, threshold = 1:10, 
#           decreasing = FALSE, GRalgo = FALSE, W = NULL, seed = NULL, 
#           as.logical = FALSE, trace = TRUE) 

#### test
# load data 

#### example
## load an MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

#### 2- data.frame function

## raw parameter
multiplot(MRIaggr.Pat1_red,param="TTP_t0",legend=FALSE,main="TTP_t0 - slice ",
             palette=rainbow(10),breaks=seq(0,10)-10^{-10})

## thresholded parameter 
data <- selectContrast(MRIaggr.Pat1_red,param=c("TTP_t0","MTT_t0","hemisphere","CSF","WM","GM"))
hypo_Th1_10 <- calcThreshold(data,param=c("TTP_t0","MTT_t0"),threshold=1:10)

multiplot(selectCoords(MRIaggr.Pat1_red),main="TTP_t0_th - slice ",
             hypo_Th1_10[,1],legend=FALSE,palette=rainbow(10),breaks=(0:10)-10^{-10})

## 1st correction
data$CSF <- as.numeric(apply(data[,c("CSF","WM","GM")],1,which.max)==1)

hypoC_Th1_10 <- calcThreshold(data,param=c("TTP_t0","MTT_t0"),threshold=1:10,
                              hemisphere="left",rm.CSF=TRUE)

multiplot(selectCoords(MRIaggr.Pat1_red),main="TTP_t0_thC - slice",
             hypoC_Th1_10[,1],legend=FALSE,palette=rainbow(10),breaks=(0:10)-10^{-10})

## 2nd correction
maskN <- c("MASK_T2_FLAIR_t2","MASK_DWI_t0")
data[,maskN] <- selectContrast(MRIaggr.Pat1_red,param=maskN)
W <- calcW(MRIaggr.Pat1_red,range=sqrt(2*1.875^2+0.001),row.norm=TRUE,upper=NULL,
                    spatial_res=c(1.875,1.875,6))
max(rowSums(W>0))

hypoCC_Th1_10 <- calcThreshold(data,param=c("TTP_t0","MTT_t0"),threshold=1:10,
                               hemisphere="left",rm.CSF=TRUE,
                               GRalgo=TRUE,seed=c("MASK_T2_FLAIR_t2","MASK_DWI_t0"),W=W)

multiplot(selectCoords(MRIaggr.Pat1_red),main="TTP_t0_thCC  - slice",
             hypoCC_Th1_10[,1],legend=FALSE,palette=rainbow(10),breaks=(0:10)-10^{-10})


#### >f calcW ####
# ?calcW
# findMethods(MRIaggr:::calcW,classes="data.frame")$data.frame
# .local <- function (object, range, method = "euclidean", 
#                     upper = NULL, format = "dgCMatrix", row.norm = FALSE, 
#                     spatial_res = rep(1, ncol(object))) 
  

#### test
data("MRIaggr.Pat1_red", package="MRIaggr")
coords <- selectCoords(MRIaggr.Pat1_red,num=1:3,hemisphere="lesion")
W <- calcW(object=coords,range=sqrt(2))
W[1:10,1:10]

W <- calcW(object=coords,range=sqrt(2),upper=TRUE)
W[1:10,1:10]


W <- calcW(object=coords,range=sqrt(2),upper=FALSE,row.norm=TRUE)
W[1:10,1:10]
table(rowSums(W))

W <- calcW(object=coords,range=sqrt(2),upper=TRUE,row.norm=TRUE)
W[1:10,1:10]
table(rowSums(W))

W <- calcW(object=coords,range=sqrt(2),spatial_res=c(1.875,1.875,6))
W[1:10,1:10]
W <- calcW(object=coords,range=6,spatial_res=c(1.875,1.875,6))
W[1:10,1:10]

W <- calcW(object=coords,range=sqrt(2),format="spam")
W[1:10,1:10]


#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

#### 1- data.frame method ####
coords <- selectCoords(MRIaggr.Pat1_red,num=1:3,hemisphere="lesion")

## full W 
W <- calcW(object=coords,range=sqrt(2))
W[1:10,1:10]
table(rowSums(W)    )

## full W normalized by row
W <- calcW(object=coords,range=sqrt(2),row.norm=TRUE)
W[1:10,1:10]
table(rowSums(W)    )

## upper W 
W <- calcW(object=coords,range=sqrt(2),upper=TRUE)
W[1:10,1:10]

#### >g EDK #### 
# pas exporte
# ?EDK
# MRIaggr:::EDK
# function (x, bandwidth, power = 2) 

#### test
bandwidth <- 2
power <- 2
x <- 10

MRIaggr:::EDK(x,bandwidth,power)
1/(sqrt(2*pi*2^bandwidth))*exp(-(x/bandwidth)^power)

#### example





#### 4- display ####

#### >a legendMRI ####
# ?legendMRI
# MRIaggr:::legendMRI
# function (breaks, palette, mar, cex, cex.main, main, quantiles,digit) 

#### test
MRIaggr:::legendMRI(1:8,palette=cm.colors(8),mar=rep(3,4),cex=1,cex.main=1,main="tot",quantiles=1:5,digit=2)

MRIaggr:::legendMRI(1:9,palette=cm.colors(9),mar=rep(3,4),cex=1,cex.main=1,main="tot",quantiles=1:5,digit=2)

#### example

#### >b multiplot ####
# ?multiplot
#  findMethods(MRIaggr:::multiplot,classes="data.frame")$data.frame
# .local <- function (object, contrast = NULL, num = NULL, 
#                     index1 = NULL, index2 = NULL, index3 = NULL, breaks = 50, 
#                     type.breaks = "range", palette = "terrain.colors", col = NULL, 
#                     pch = 15, cex = 1, col.NA = "lightyellow", pch.NA = 8, 
#                     xlim = NULL, ylim = NULL, axes = TRUE, window = FALSE, 
#                     legend = TRUE, mfrow = NULL, mar = rep(1.5, 4), mgp = c(2,0.5, 0), pty = NULL, asp = 1, bg = "lightblue", xlab = "", 
#                     ylab = "", main = "slice ", num.main = TRUE, cex.main = 0.9, 
#                     quantiles.legend = TRUE, digit.legend = 3, cex.legend = 1.5, 
#                     mar.legend = c(2, 7, 2, 2), main.legend = names(contrast), 
#                     filename = "multiplot", width = 1000, height = 700, 
#                     path = NULL, unit = "px", res = NA) 

#### test
data <- selectContrast(MRIaggr.Pat1_red,param=c("DWI_t0","TTP_t0","MTT_t0","MASK_T2_FLAIR_t2"),
                        hemisphere="lesion",coords=TRUE)
glm.1 <- glm(MASK_T2_FLAIR_t2 ~ DWI_t0 + TTP_t0 + MTT_t0,data=data,family=binomial(link="logit"))

multiplot(object=data[,c("i","j","k")],num=2:3,
             contrast=predict(glm.1,type="response"),window=FALSE,
             quantiles.legend=FALSE,
             index1=list(coords=data[data$MASK_T2_FLAIR_t2,c("i","j","k")],outline=TRUE)
)

multiplot(object=data[,c("i","j","k")],num=5,
             contrast=predict(glm.1,type="response"),window=FALSE,
             type.breaks="range_center",col=rgb(predict(glm.1,type="response"),0,0),
             index1=list(coords=data[data$MASK_T2_FLAIR_t2,c("i","j","k")],outline=TRUE)
)

data.test <- rbind(data[,c("i","j","k")],cbind(data[,c("i","j")],k=data$k+6))
multiplot(object=data.test,
             contrast=c(predict(glm.1,type="response"),predict(glm.1,type="response")))
multiplot(object=data.test,window="png",path="E:/",
             contrast=c(predict(glm.1,type="response"),predict(glm.1,type="response")))

multiplot(object=data[,c("i","j","k")],num=1,
             contrast=c(rep(NA,100),predict(glm.1,type="response")[-(1:100)]),
             type.breaks="range_center",col=rgb(predict(glm.1,type="response"),0,0)
)

multiplot(object=data[,c("i","j","k")],num=1,
             contrast=c(rep(NA,100),predict(glm.1,type="response")[-(1:100)]),window=FALSE,
             type.breaks="range_center",col=rgb(predict(glm.1,type="response"),0,0),
             index1=list(coords=data[data$MASK_T2_FLAIR_t2,c("i","j","k")],outline=TRUE)
)


multiplot(object=data[,c("i","j","k")],num=1,
             contrast=c(rep(NA,100),predict(glm.1,type="response")[-(1:100)]),
             type.breaks="range_center")

multiplot(object=data[,c("i","j","k")],num=1,
             contrast=c(rep(NA,100),predict(glm.1,type="response")[-(1:100)]),
             type.breaks="range_center",col.NA=NA)


#### example
## simulate 
n <- 10
Y <- rnorm(n^2)

## display
multiplot(object=data.frame(expand.grid(1:n,1:n),1),
             contrast=Y,window=FALSE)

## load an MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

## select data
data <- selectContrast(MRIaggr.Pat1_red,param=c("DWI_t0","TTP_t0","MTT_t0","MASK_T2_FLAIR_t2"),
                    hemisphere="lesion",coords=TRUE)

## fit model
glm.1 <- glm(MASK_T2_FLAIR_t2 ~ DWI_t0 + TTP_t0 + MTT_t0,data=data,family=binomial(link="logit"))

## display fitted values
multiplot(object=data[,c("i","j","k")],
             contrast=predict(glm.1,type="response"),window=FALSE)

## display residuals
multiplot(object=data[,c("i","j","k")],num=3,
             contrast=predict(glm.1,type="response"),window=FALSE,
             index1=list(coords=data[data$MASK_T2_FLAIR_t2,c("i","j","k")],outline=TRUE)
)

#### >c outline ####
# ?outline
# outline
# function (n = 50, sequential = FALSE, min_dist = 1, col = c("blue","red", "grey"), pch = 20, cex = c(0.75, 1, 0.75))

#### test

#### example
## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")
num <- 3

## display 1
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",              
             num=num,legend=FALSE,             
             window=FALSE)

## outline on display 1
res <- outline(sequential=TRUE,min_dist=3)

## display the results
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",              
             num=num,legend=FALSE,
             index1=data.frame(k=num,res$edge[,c("i","j")]),
             index2=data.frame(k=num,res$surface[,c("i","j")]),
             window=FALSE)

carto <- selectContrast(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","index"),num=num,coords=TRUE)
carto <- merge(carto,cbind(res$surface,outline=TRUE),all = TRUE)
carto[is.na(carto$outline),"outline"] <- FALSE
head(carto)

## display the results next to MASK_T2_FLAIR_t2
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",              
             num=num,legend=FALSE,
             index1=carto[carto$MASK_T2_FLAIR_t2,c("i","j","k")],
             index2=carto[carto$outline,c("i","j","k")],
             window=FALSE)


#### >d plotMRI ####
# ?plotMRI
# MRIaggr:::plotMRI
# function (contrast, coords, breaks, palette, col, plot.image, asp, 
#           xlim, ylim, pch, cex, axes, col.NA, pch.NA, xlab, ylab, main, 
#           cex.main) 

#### test

#### example

set.seed(10)
n <- 10
Y <- rnorm(n^2)
res <- MRIaggr:::plotMRI(contrast=as.matrix(Y), coords=expand.grid(1:n,1:n), 
                         breaks=seq(-10,10), palette=terrain.colors(20), col=NULL, 
                         asp=1, xlim=c(0,n+1), ylim=c(0,n+1), 
                         pch=NULL, cex=1, axes=TRUE, col.NA="yellow", pch.NA=20, 
                         xlab="", ylab="", main="", cex.main=1) 

#### >e pointsOutline ####
# ?pointsOutline
# pointsOutline
# function (coords, array = NULL, filter = "2D_N4") 

#### test
edge <- pointsOutline(expand.grid(1:5,1:5))

plot(expand.grid(1:5,1:5))
points(edge,col="red",pch=20)

#### example


#### 5- init ####

#### >a initCol  ####
# ?initCol
# MRIaggr:::initCol
# function (contrast, param = NULL, pch, col, palette, breaks, legend, type.breaks, 
#           method = "initCol")

#### test
MRIaggr:::initCol(contrast=cbind(1:10), coords=1:10, pch=NULL, param = "TMAX", col=NULL, palette="terrain.colors", breaks=50, legend=FALSE, type.breaks="quantile", 
          method = "initCol")
#### example

#### >b initDisplayWindow  ####
# ?initDisplayWindow
# MRIaggr:::initDisplayWindow
# function (window, filename, path, width, height, scale, res, 
#           mfrow, bg, pty, mar, mgp) 

#### test
MRIaggr:::initDisplayWindow(window="png", filename="test", path=NULL, width=400, height=400, scale=1, res=NA, 
                            mfrow=c(1,1), bg="red", pty="m", mar=rep(2,4), mgp=c(2,1,0)) 
dev.off()

#### example

#### >c initFilter  ####
# ?initFilter
# initFilter
# function(filter,method="calcFilter")

#### test
#### examples
initFilter("2D_G3",method="calcFilter")$filter
initFilter("2D_G5",method="calcFilter")$filter
initFilter("3D_G3",method="calcFilter")$filter

initFilter("2D_M9",method="calcFilter")$filter
initFilter("3D_M3",method="calcFilter")$filter

initFilter("2D_Sx",method="calcFilter")$filter
initFilter("2D_Sy",method="calcFilter")$filter
initFilter("3D_Sx",method="calcFilter")$filter
initFilter("3D_Sz",method="calcFilter")$filter

initFilter("2D_I3",method="calcFilter")$filter
initFilter("2D_I5",method="calcFilter")$filter
initFilter("2D_I7",method="calcFilter")$filter
initFilter("2D_I9",method="calcFilter")$filter
initFilter("3D_I5",method="calcFilter")$filter


#### >d initIndex  ####
# ?initIndex
# MRIaggr:::initIndex
# function (object, index, num, hemisphere = "both", as.logical, 
# indexNum = NULL, cex.default, pch.default, col.default, filter_default, 
# method = "initIndex") 

#### test
data("MRIaggr.Pat1_red", package="MRIaggr")
res <- MRIaggr:::initIndex(object=MRIaggr.Pat1_red, index="MASK_DWI_t0", 
          num=NULL, hemisphere = "both", as.logical=TRUE, 
          indexNum = NULL, cex.default=1, pch.default=15, col.default="red", filter_default="2D_4N", 
          method = "initIndex") 

dfex <- data.frame(expand.grid(1:10, 1:10),1:5)
names(dfex) <- c("i","j","k")
res <- MRIaggr:::initIndex(object=data.frame(), index=dfex, 
                 num=NULL, hemisphere = NULL, as.logical=FALSE, 
                 indexNum = NULL, cex.default=1, pch.default=15, col.default="red", filter_default="2D_4N", 
                 method = "initIndex") 

#### example

####  >e initNeighborhood  ####
# ?initNeighborhood
# initNeighborhood
# function (Neighborhood_name, method = "initNeighborhood") 

#### test
#### examples
# 2D neighborhood
initNeighborhood("2D_N4",method="calcFilter") # rock neighborhood
initNeighborhood("2D_N8",method="calcFilter") # queen neighborhood

# 3D neighborhood
initNeighborhood("3D_N6",method="calcFilter") # rock neighborhood
initNeighborhood("3D_N10",method="calcFilter") # queen neighborhood
initNeighborhood("3D_N18",method="calcFilter") # queen neighborhood
initNeighborhood("3D_N26",method="calcFilter") # queen neighborhood

#### >f initWindow  ####
# ?initWindow
# MRIaggr:::initWindow
# function (window, filename, path, width, height, unit, res, n.plot, 
# mfrow, xlim, ylim, method = "initWindow") 

#### test
MRIaggr:::initWindow(window="png", filename="test", path=NULL, width=400, height=400, unit="px", res=NA, n.plot=3, 
          mfrow=NULL, xlim=NULL, ylim=NULL, method = "initWindow") 

#### example

