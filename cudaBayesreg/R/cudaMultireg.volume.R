## Run cudaMultireg.slice on a data volume

cudaMultireg.volume <-
function(fbase=NULL, R=2000, keep=5, nu.e=3, zprior=FALSE, rng=0, rg=c(NULL,NULL),
  swap=FALSE, savedir=tempdir())
{
  if(!is.null(fbase)) {
    fsl.filtered <- system.file(paste("extdata/", fbase,
      "_filtered_func_data.nii.gz", sep = ""), package = "cudaBayesregData")
    fsl.mask <- system.file(paste("extdata/", fbase, "_mask.nii.gz", 
      sep = ""), package = "cudaBayesregData")
    ## Design matrix (FSL-like) 
    fsl.design <- system.file(paste("extdata/", fbase, "_design.mat", 
      sep = ""), package = "cudaBayesregData")
  } else {
    fsl.filtered <- "filtered_func_data.nii.gz"
    fsl.mask <- "mask.nii.gz"
    fsl.design <- "design.mat"
  }
  img.nifti <- readNIfTI(fsl.filtered)
  img <- img.nifti@.Data
	d <- dim(img)
	if(is.null(rg)) {
		nslices <- d[3]
		first <- 1
		last <- nslices
	}
	else {
		first <- rg[1]
		last <- rg[2]
	}
	for (sl in (first:last)) {
		cat("\n*** multilevel slice n.", sl,"\n")
		slicedata <- read.fmrislice(fbase=fbase, slice=sl, swap=swap)
		ymaskdata <- premask(slicedata)
		fsave <- paste(savedir,"/",fbase,"_s",sl,"_nu",nu.e,".sav",sep="")
		out <- cudaMultireg.slice(slicedata, ymaskdata, R=R, keep=keep, nu.e=nu.e,
		  fsave=fsave, zprior=zprior, rng=rng)
	}
}

