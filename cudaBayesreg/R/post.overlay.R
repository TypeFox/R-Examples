#
# rendering zstat overlays
#
post.overlay <-
function(fbase=NULL, vreg=2, nu.e=3, rg=c(NULL,NULL), view="axial", savedir=tempdir())
{
    if(!is.null(fbase)) {
        fsl.filtered <- system.file(paste("extdata/", fbase,
          "_filtered_func_data.nii.gz", sep = ""), package = "cudaBayesregData")
    }
    else {
        fsl.filtered <- "filtered_func_data.nii.gz"
    }
#####################
#  fsl.filtered <- system.file(paste("extdata/", fbase,
#    "_filtered_func_data.nii.gz", sep = ""), package = "cudaBayesregData")
#####################
	a <- readNIfTI(fsl.filtered)@.Data[,,,1]
	cat(paste("loaded", fsl.filtered ,"\n"))
	da <- dim(a)
	zstatfname <- paste(savedir,"/",fbase,"_zstat",vreg,"_nu",nu.e,sep="")
	b <- readNIfTI(zstatfname)@.Data
	cat(paste("loaded", zstatfname ,"\n"))
	db <- dim(b)
	stopifnot(da == db)
	a1 <- da[1]; a2 <- da[2]; a3 <- da[3];
	##--------
	normalize <- function(xv) { # normalize data 
	 	r <- range(xv); rd <- r[2]-r[1];
	 	if(rd) {xv <- (xv-r[1])/rd};
	 	invisible(xv)
	}
	a <- normalize(a)
	b <- normalize(b)
	# a <- 0.8 * normalize(a)
	# b <- 0.2 * normalize(b)
	##--------
	if(is.null(rg)) 
		rg <- switch(EXPR = view, axial=c(1,a3), coronal=c(1,a2), sagittal=c(1,a1))
	nfig <- rg[2]-rg[1] + 1
	nc <- ceiling(sqrt(nfig))
	nr <- ceiling(nfig/nc)
	zlim <- range(b)
	b [ b == 0] <- NA 
	par(mfrow=c(nr,nc),par(mar=c(0,0,0,0)+0.1))
	for(i in rg[1]:rg[2]) {
		ai <- switch(EXPR = view,
			axial=a[,,i], coronal=a[,i,], sagittal=a[i,,])
		bi <- switch(EXPR = view,
			axial=b[,,i], coronal=b[,i,], sagittal=b[i,,])
		image(ai, col=gray((0:255)/256), xaxt="n", yaxt="n")
		image(bi, col=heat.colors(256), xaxt="n", yaxt="n", add=TRUE, zlim=zlim)
	}
}






