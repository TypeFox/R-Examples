#
# Mask out slice times series and keep indices
#
premask <-
function (slicedata) 
{
    slice <- slicedata$slice
    niislicets <- slicedata$niislicets
    mask <- slicedata$mask
    kin <- which(mask == 1, arr.ind = T) # indices of pixels in mask 
    if (!length(kin)) {
        cat("\n slice", slice, ":\tempty slice mask - nothing to do\n")
        return()
    }
    d <- dim(kin)
    ym <- NULL
    for (i in 1:d[1]) {
        yx <- niislicets[kin[i, 1], kin[i, 2], ]
        if (sd(yx)) {  # do not include null time series even if mask is 1 
            ym <- cbind(ym, yx)
        }
        else { # remove form mask 
            mask[kin[i, 1], kin[i, 2]] <- 0
        }
    }
    kin <- which(mask == 1, arr.ind = T)  # update indices of pixels in mask
		###
    stdf <- function(y) { return((y - mean(y))/sd(y)); }
    yn <- apply(ym, 2, stdf)
    nobs <- slicedata$nobs
    stopifnot(nobs == nrow(yn))
    nreg <- ncol(yn)
    invisible(list(yn = yn, kin = kin, nreg = nreg))
}
