AvBasis.wst2D <-
function (wst2D, ...) 
{
    filter <- wst2D$filter
    amdim <- dim(wst2D$wst2D)
    im <- matrix(0, nrow = amdim[2]/2, ncol = amdim[2]/2)
    ans <- .C("SAvBasis", am = as.double(wst2D$wst2D), d1 = as.integer(amdim[1]), 
        d12 = as.integer(amdim[1] * amdim[2]), TheSmooth = as.double(im), 
        levj = as.integer(amdim[1]), H = as.double(filter$H), 
        LengthH = as.integer(length(filter$H)), error = as.integer(0), 
        PACKAGE = "LS2W")
    if (ans$error != 0) 
        stop(paste("Error code was ", ans$error))
    matrix(ans$TheSmooth, nrow = amdim[2]/2)
}
