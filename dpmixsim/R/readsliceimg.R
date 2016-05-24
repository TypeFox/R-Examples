## MRI slice data as prefiltered by FSL
## 	niislice: slice image
##		mask: slice mask

readsliceimg <-
function (fbase="t1_pn3_rf0", swap=FALSE) 
{
    mri.filtered <- system.file(paste("extdata/", fbase, "_slice_0092.nii.gz", 
        sep = ""), package = "dpmixsim")
    mri.mask <- system.file(paste("extdata/", fbase, "_slice_0092_mask.nii.gz", 
        sep = ""), package = "dpmixsim")
    img.nifti <- readNIfTI(mri.filtered)
    mask.nifti <- readNIfTI(mri.mask)
    img <- img.nifti@.Data
    mask <- mask.nifti@.Data
    nr <- nrow(img);
		nc <- ncol(img)
    nrm <- nrow(mask);
		ncm <- ncol(mask)
		stopifnot(nr == nrm)
		stopifnot(nc == ncm)
    if (swap) { 
        niislice <- img[nr:1,]
        mask <- mask[nr:1,]
    }
    else {
        niislice <- img
        mask <- mask
    }
    invisible(list(fbase=fbase, niislice=niislice, mask=mask, nrow=nr, ncol=nc, swap=swap))
}
