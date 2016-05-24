#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Fichier test 2.1 : Test des fonctions associees a Carto3D
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require(MRIaggr)

options(error=function() traceback(2)) 
options(max.print=10000)

#### 0 chargement ####
path.nifti_files <- system.file(file.path("nifti"),package = "MRIaggr")

nifti.Pat1_TTP_t0 <- readMRI(file=paste(path.nifti_files,"TTP_t0",sep="/"),
                             format="nifti")


Carto3D.Pat1_TTP_t0.new <- new("Carto3D",
                               identifier = "Pat1",
                               parameter	= "TTP_t0",
                               contrast	= array2df(nifti.Pat1_TTP_t0[,,,1],name_newparam="TTP_t0"),
                               voxelDim = data.frame(i=128,j=128,k=6),
                               default_value	= as.character(nifti.Pat1_TTP_t0[1,1,1,1])
)

Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,
                                    identifier="Pat1",param="TTP_t0",default_value=NULL,
                                    pos_default_value=cbind(1:10,1,1))

#### 1a selectContrast ####

#### description
# ?selectContrast
# findMethods(selectContrast,classes="Carto3D")$Carto3D
# local <- function (object, num = NULL, na.rm = FALSE, coords = TRUE, 
#                    format = "any") 

#### test
res <- selectContrast(Carto3D.Pat1_TTP_t0)
res <- selectContrast(Carto3D.Pat1_TTP_t0,num=1)
res <- selectContrast(Carto3D.Pat1_TTP_t0,num=1:3)

res <- selectContrast(Carto3D.Pat1_TTP_t0,coords=FALSE)
res <- selectContrast(Carto3D.Pat1_TTP_t0,num=1,coords=FALSE)
res <- selectContrast(Carto3D.Pat1_TTP_t0,format="matrix",coords=FALSE)
res <- selectContrast(Carto3D.Pat1_TTP_t0,format="matrix",coords=FALSE,na.rm=TRUE)

#### example 
## load NIFTI files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## select all observations
carto1 <- selectContrast(Carto3D.Pat1_TTP_t0)
dim(carto1)

## select observations from slices 1 to 3 and return the result into a data.frame
carto2 <- selectContrast(Carto3D.Pat1_TTP_t0,num=1:3,coords=FALSE,format="data.frame")
dim(carto2)

## select observations from slices 1 to 3 and return the result into a vector
carto3 <- selectContrast(Carto3D.Pat1_TTP_t0,num=1:3,coords=FALSE)
length(carto3)

#### 1b selectCoords ####

#### description
# ?selectCoords
# findMethods(selectCoords,classes="Carto3D")$Carto3D
# .local <- function (object, coords = c("i", "j", "k"), num = NULL, 
#                     format = "any") 


#### test
head(selectCoords(Carto3D.Pat1_TTP_t0))

res <- selectCoords(Carto3D.Pat1_TTP_t0,num=1,coords="i")
res <- selectCoords(Carto3D.Pat1_TTP_t0,num=1,coords="i",format="data.frame")
res <- selectCoords(Carto3D.Pat1_TTP_t0,num=1)

#### example
## load NIFTI files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## selection all coordinates
coords1 <- selectCoords(Carto3D.Pat1_TTP_t0)
dim(coords1)

## selection coordinates i and j from slices 1 to 3
coords2 <- selectCoords(Carto3D.Pat1_TTP_t0,num=1:3,coords=c("i","j"))
dim(coords2)


#### 1c selectDefault_value ####

#### description
# ?selectDefault_value
# findMethods(selectDefault_value,classes="Carto3D")$Carto3D
# .local <- function (object) 

#### test

#### example
## load NIFTI files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## selection
selectDefault_value(Carto3D.Pat1_TTP_t0)

#### 1d selectVoxelDim ####

#### description
# ?selectVoxelDim
# findMethods(selectVoxelDim,classes="Carto3D")$Carto3D
# .local <- function (object) 
  
#### test

#### example
## load nifti files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## selection
selectVoxelDim(Carto3D.Pat1_TTP_t0)

#### 1e selectIdentifier ####

#### description
# ?selectIdentifier
# findMethods(selectIdentifier,classes="Carto3D")$Carto3D
# .local <- function (object) 
  
#### test

#### example
## load NIFTI files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## selection
selectIdentifier(Carto3D.Pat1_TTP_t0)

#### 1f selectParameter ####

#### description
# ?selectParameter
# findMethods(selectParameter,classes="Carto3D")$Carto3D
# .local <- function (object) 

#### test

#### example
## load NIFTI files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## selection
selectParameter(Carto3D.Pat1_TTP_t0)

#### 3a multiplot ####

#### description
# ?multiplot
# findMethods(multiplot,classes="Carto3D")$Carto3D
# .local <- function (object, num = NULL, breaks = 50, type.breaks = "range", 
#                     palette = "terrain.colors", col = NULL, pch = 15, cex = 1, 
#                     col.NA = "lightyellow", pch.NA = 8, xlim = NULL, ylim = NULL, 
#                     axes = TRUE, window = FALSE, legend = FALSE, mfrow = NULL, 
#                     mar = rep(1.5, 4), mgp = c(2, 0.5, 0), pty = NULL, asp = 1, 
#                     bg = "lightblue", xlab = "", ylab = "", main = NULL, 
#                     num.main = TRUE, cex.main = 1.5, quantiles.legend = TRUE, 
#                     digit.legend = 3, cex.legend = 1.5, mar.legend = c(1,4, 2, 1), filename = "multiplot", width = 1000, 
#                     height = 700, path = NULL, unit = "px", res = NA) 

#### test
Carto3D.Pat1_TTP_t0@contrast[,4] <- 1
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=F,main="",num.main=FALSE,
             palette="gray.colors",breaks=seq(0,100))

multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=F,main="",num.main=FALSE,
             palette="red")

Carto3D.Pat1_TTP_t0@contrast[,4] <- c(1,2)

multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=F,main="",num.main=FALSE,
             palette=c("red","blue"))

Carto3D.Pat1_TTP_t0@contrast[,4] <- c(1,2,3)

multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=F,main="",num.main=FALSE,
             palette=c("red","blue","green"))

#### example
## load NIFTI files and convert them to carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## display
multiplot(Carto3D.Pat1_TTP_t0)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE,legend=FALSE)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE,legend=FALSE,main="",num.main=FALSE)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE,main="",num.main=FALSE,
             palette="gray.colors",breaks=seq(0,100))

#### 5a initNum ####

#### description
# ?MRIaggr:::initNum
# findMethods(MRIaggr:::initNum,classes="Carto3D")$Carto3D
# .local <- function (object, num, test = TRUE, init = TRUE) 

#### test
MRIaggr:::initNum(Carto3D.Pat1_TTP_t0,num=1:3)
MRIaggr:::initNum(Carto3D.Pat1_TTP_t0,num=NULL)
MRIaggr:::initNum(Carto3D.Pat1_TTP_t0,num=NULL,init=FALSE)
MRIaggr:::initNum(Carto3D.Pat1_TTP_t0,num=-1,test=FALSE)

#### example

