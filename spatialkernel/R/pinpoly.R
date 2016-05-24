## -1 outside; 0 on boundary; 1 inside
## number of poly points more than 3000--error
## return values change to -1-error, 0-outside, 1-boundary, 2-inside
pinpoly<-function(poly, pts)
{
    if(nrow(poly)>3000) {
        cat("\nBoundary polygon vertices number exceeds 3000.\n")
        return(-1)
    }
    if(is.matrix(pts)){
        ans<-.Fortran("psnpoly", as.double(pts[,1]), as.double(pts[,2]),
                      as.integer(nrow(pts)), as.double(poly[,1]), as.double(poly[,2]),
		      as.integer(nrow(poly)), 
                      inout=integer(nrow(pts)), PACKAGE="spatialkernel")$inout
    }else{
        ans<-.Fortran("pnpoly", as.double(pts[1]), as.double(pts[2]),
                      as.double(poly[,1]), as.double(poly[,2]), as.integer(nrow(poly)), 
                      inout=as.integer(0), PACKAGE="spatialkernel")$inout
    }
    ans + 1
}
