path <- system.file("nifti",package = "MRIaggr")
require(MRIaggr)
data("MRIaggr.Pat1_red", package="MRIaggr")
system.time(
  res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",
                        trace=TRUE,update.object=TRUE,overwrite=TRUE)
)
multiplot(MRIaggr.Pat1_red,param="DWI_t0",midplane=TRUE)
test.package <- requireNamespace("EBImage",quietly=TRUE)
if(test.package){
## load data and build MRIaggr
MRIparameters <- c("T2_GRE_t0","MASK_DWI_t0","TTP_t1")
list.nifti <- list()
for(iter_param in 1:length(MRIparameters)){
  iter_file <- file.path(path,MRIparameters[iter_param])
  list.nifti[[iter_param]] <- readMRI(file=iter_file,format="nifti")
}

voxelDim <- data.frame(i=1.875,j=1.875,k=6,unit="mm")
MRIaggr.Pat1 <- constMRIaggr(list.nifti,identifier="Pat1",
                             param=MRIparameters,voxelDim=voxelDim)

## image rotation
data.df <- selectContrast(MRIaggr.Pat1,param=MRIparameters)
coords.df <- selectCoords(MRIaggr.Pat1)
data.array <- df2array(contrast=data.df,coords=coords.df)

data.array_flip <- EBImage::rotate(data.array[["contrast"]][["T2_GRE_t0"]],30)[25:(24+78),13:(12+100),1:3]
# image(data.array_flip[,,2])
data.df_flip <- array2df(data.array_flip)
allocContrast(MRIaggr.Pat1,param="T2_GRE_t0_flip")  <- data.df_flip

## create the cerebral mask
res <- calcBrainMask(MRIaggr.Pat1,param="T2_GRE_t0_flip",type="kmeans",
                     kmeans.n_groups=2:4,
                     update.object=TRUE,overwrite=TRUE)

## smooth the cerebral mask
res <- calcSmoothMask(MRIaggr.Pat1,update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1,param="mask",legend=FALSE)

## construct the reduced object
MRIaggr.Pat1_red <- constReduceMRIaggr(MRIaggr.Pat1,mask="mask")

## both
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0_flip",
                      mask=c("MASK_DWI_t0"),as.logical=TRUE,
                      i_test=seq(-10,10,1),angle_test=seq(-40,-20,1),
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0_flip",num=1,legend=FALSE,
          midplane=TRUE,main="original coordinates - slice ")

deg2rad <- 2*pi/360
0.52/deg2rad

## grid Search only
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0_flip",
                      mask=c("MASK_DWI_t0"),as.logical=TRUE,
                      NelderMead=FALSE,
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0_flip",num=1,legend=FALSE,
          midplane=TRUE,main="original coordinates - slice ")


## NelderMead only
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0_flip",
                      mask=c("MASK_DWI_t0"),as.logical=TRUE,
                      gridSearch=FALSE,
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0_flip",num=1,legend=FALSE,          
          midplane=TRUE,main="original coordinates - slice ")

## no rotation
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",
                      mask=c("MASK_DWI_t0"),as.logical=TRUE,
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)

multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,legend=FALSE,
          midplane=TRUE,main="original coordinates - slice ")



demo("preProcessing")


## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## Not run: 
system.time(
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)
)
## End(Not run)


## display the mid-sagittal plan
multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0",num=3,legend=FALSE,
          midplane=TRUE,main="original coordinates - slice ")

## display
multiplot(selectContrast(MRIaggr.Pat1_red,param=c("i_hemisphere","j_hemisphere","k")),
          contrast=selectContrast(MRIaggr.Pat1_red,param="T2_GRE_t0"),num=3,
          index1=cbind(0,seq(-50,50),3),main="new coordinates - slice ",legend=FALSE)


## compute the mid-sagittal plan and mark lesion
## Not run: 
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",
                      mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)
}