#
# read segmented images to build Z array
# Z prior as csf/gry/wht segmented regions 
#	Extensions for segmented masks as in *_csf.nii.gz, etc.
#
read.Zsegslice <-
function(slicedata, ymaskdata)
{
	fbase <- slicedata$fbase
	slice <- slicedata$slice
  if(!is.null(fbase)) {
    fsl.csf <- system.file(paste("extdata/",fbase,"_csf.nii.gz",sep=""),
		  package = "cudaBayesregData")
    fsl.gry <- system.file(paste("extdata/",fbase,"_gry.nii.gz",sep=""),
		  package = "cudaBayesregData")
    fsl.wht <- system.file(paste("extdata/",fbase,"_wht.nii.gz",sep=""),
	  	package = "cudaBayesregData")
  }
  else {
    fsl.csf <- "csf.nii.gz"
    fsl.gry <- "gry.nii.gz"
    fsl.wht <- "wht.nii.gz"
  }
  stopifnot(file.exists(fsl.csf), file.exists(fsl.gry), file.exists(fsl.wht))
  csfm.nifti <- readNIfTI(fsl.csf)
  csfm <- csfm.nifti@.Data
	grym.nifti <- readNIfTI(fsl.gry)
  grym <- grym.nifti@.Data
	whtm.nifti <- readNIfTI(fsl.wht)
  whtm <- whtm.nifti@.Data
	nx <- nrow(csfm)
	swap <- slicedata$swap
  if(swap) { # read Z-mask slice data consistent with slicedata
		csf.sl <- csfm[nx:1,,slice]
		gry.sl <- grym[nx:1,,slice]
		wht.sl <- whtm[nx:1,,slice]
	} else {
		csf.sl <- csfm[,,slice]
		gry.sl <- grym[,,slice]
		wht.sl <- whtm[,,slice]
	}
	#----------------
	kin <- ymaskdata$kin
	nreg <- ymaskdata$nreg
	d <- dim(kin) 
	zmat <- matrix(0,nrow=nreg,ncol=3)
	for(i in 1:d[1]) {
		c0 <- csf.sl[kin[i,1],kin[i,2]]
		zmat[i,1] <- c0
		c0 <- gry.sl[kin[i,1],kin[i,2]]
		zmat[i,2] <- c0
		c0 <- wht.sl[kin[i,1],kin[i,2]]
		zmat[i,3] <- c0
	}
	demean <- function(x) { x-mean(x); }
	zmat <- apply(zmat,2,demean)
	Z <- cbind(c(rep(1,nreg)), zmat)
	invisible(Z)
}
