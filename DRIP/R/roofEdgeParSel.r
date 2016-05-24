# This is R source code for function 'roofEdgeParSel' in
# the R package 'image'.
# Creator: Yicheng Kang
# Date: Sep 13, 2015

roofEdgeParSel = function(image, bandwidth, thresh, nboot, edge1, blur=FALSE) {
    if (!is.matrix(image))
        stop("image data must be a matrix")
    n1 = as.integer(dim(image)[1])
    n2 = as.integer(dim(image)[2])
    bandwidth = as.integer(bandwidth)
    nboot = as.integer(nboot)
    nband = length(bandwidth)
    if (n1 != n2)
        stop("image data must be a square matrix")
    if (!is.numeric(bandwidth))
        stop("bandwidth must be numeric")
    if (min(bandwidth) < 1)
        stop("every candidate bandwidth must be a positive integer")
    if (n1 + 2 * max(bandwidth) + 2 > 600)
        stop("some bandwidth is too large or the resolution
of the image is too high.")
    if (!is.numeric(thresh))
        stop("threshold must be numeric")
    if ((!is.numeric(nboot)) | (length(nboot) > 1) |
        (as.integer(nboot) < 1))
        stop('nboot must be a positive integer.')
    if(!is.matrix(edge1) | ncol(edge1) != nrow(edge1))
        stop("edge1 must be a square matrix")
    if(!all(edge1==0 | edge1==1))
        stop("edge1's must be either 0 or 1.")
    if(ncol(edge1) != n1)
        stop("edge1 and image are not of the same size.")
    n = dim(image)[1]
    z = matrix(as.double(image), ncol = n)
    edge1 = matrix(as.integer(edge1), ncol=n)
    u = as.double(thresh)
    nthresh = length(u)
    out.mat = array(as.double(0), c(nband, nthresh))
    if (blur == FALSE) {
        out = .Fortran("roofEdgeParSel_denoise", n=as.integer(n-1), obsImg=z, nband=nband,
            bandwidth=bandwidth, nthresh=nthresh, thresh=u, nboot=nboot, edge1=edge1,
            dKQ=out.mat)
    }
    else {
        out = .Fortran("roofEdgeParSel_deblur", n=as.integer(n-1), obsImg=z, nband=nband,
            bandwidth=bandwidth, nthresh=nthresh, thresh=u, nboot=nboot, edge1=edge1,
            dKQ=out.mat)
    }
    out.mat = out$dKQ
    rownames(out.mat) = paste('bandwidth=', bandwidth, sep='')
    colnames(out.mat) = paste('thresh=', thresh, sep='')
    out.mat.min = min(out.mat)
    for(iband in 1:nband){
    	for(ithresh in 1:nthresh){
            if(out.mat[iband, ithresh] == out.mat.min){
                band_sel = bandwidth[iband]
                thresh_sel = thresh[ithresh]
            }
    	}
    }
    paste('The selected bandwidth is', band_sel)
    paste('The selected threshold is', thresh_sel)
    return(list(output_matrix=out.mat, selected_bandwidth=band_sel, selected_threshold=thresh_sel))
}
