## fMRI data as prefiltered by FSL :  read all time series in slice  
##   niislicets: slice time-series
##    mask: slice mask
##   X: design matrix
##    nvar: n. of regression vars
##    nobs: n. of observations

read.design <-
function(fsldesign="design.mat")
{
    g1   <- paste("grep NumWaves", fsldesign, " | awk '{ print $2 }'")
    nevs <-as.integer(system(g1, intern=TRUE))
    stopifnot(!is.na(nevs))
    g2   <- paste("grep NumPoints", fsldesign," | awk '{ print $2 }'")
    npts <- as.integer(system(g2, intern=TRUE))
    stopifnot(!is.na(npts))
    ## cat("npts=",npts,"\n")
    ## cat("nevs=",nevs,"\n")
    dsgn <- numeric(nevs*npts)
    model <- .C("readmodel",
      fsldesign=as.character(fsldesign),
      nevs=as.integer(nevs),
      npts=as.integer(npts),
      dsgn=as.double(dsgn) )
    invisible(model)
}

read.fmrislice <-
function (fbase=NULL, slice=NULL, swap=FALSE) 
{
    if(!is.null(fbase)) {
        fsl.filtered <- system.file(paste("extdata/", fbase,
          "_filtered_func_data.nii.gz", sep = ""), package = "cudaBayesregData")
        fsl.mask <- system.file(paste("extdata/", fbase, "_mask.nii.gz", 
          sep = ""), package = "cudaBayesregData")
        ## Design matrix (FSL-like) 
        fsl.design <- system.file(paste("extdata/", fbase, "_design.mat", 
          sep = ""), package = "cudaBayesregData")
    }
    else {
        fsl.filtered <- "filtered_func_data.nii.gz"
        fsl.mask <- "mask.nii.gz"
        fsl.design <- "design.mat"
    }
    stopifnot(file.exists(fsl.design))
    img.nifti <- readNIfTI(fsl.filtered)
    mask.nifti <- readNIfTI(fsl.mask)
    img <- img.nifti@.Data
    mask <- mask.nifti@.Data
    dimg <- img.nifti@dim
    if(is.null(slice)) ## default slice n.
       slice <- floor(dimg[3]/2)
    cat("Using slice n.", slice,"\n")
    X <- nrow(img)
    Xm <- nrow(mask)
    if (swap) { ## swap=T to be represented as in fslview
        niislicets <- img[X:1, , slice, ]
        mask <- mask[Xm:1, , slice]
    }
    else {
        niislicets <- img[, , slice, ]
        mask <- mask[, , slice]
    }

    model <-read.design(fsldesign=fsl.design)
    dsgn <- matrix(model$dsgn, ncol=model$nevs, byrow=T)
    nobs <- nrow(dsgn)
    X0 <- as.matrix(dsgn)
    X <- cbind(rep(1, nobs), X0) ## with intercept
    nvar <- ncol(X)
    invisible(list(fbase=fbase, slice = slice, niislicets = niislicets, mask = mask, 
        X = X, nvar = nvar, nobs = nobs, swap = swap))
}



