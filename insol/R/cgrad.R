cgrad <-
function(dem, dlx=0, dly=dlx, cArea=FALSE){
    if (nargs() < 1) {
        cat("USAGE: cgrad(dem, dx, dly=dlx, cArea=FALSE) \n")
        return()
    }
	if (class(dem)=="RasterLayer") {
		dlx = res(dem)[1]
		dly = res(dem)[2]
		dem=raster::as.matrix(dem)
	}
	if (dlx == 0){
				cat("Input data is not a RasterLayer, then I need the DEM resolution dlx \n")
		return()
		}
	mm=as.matrix(dem)
	rows=nrow(mm)
	cols=ncol(mm)
	cellgr=array(dim=c(rows,cols,3))
	md=mm[-rows,-1]
	mr=mm[-1,-cols]
	mrd=mm[-1,-1]
	cellgr[-rows,-cols,2]=.5*dlx*(mm[-rows,-cols]+md-mr-mrd)
	cellgr[-rows,-cols,1]=.5*dly*(mm[-rows,-cols]-md+mr-mrd)
	cellgr[-rows,-cols,3]=dlx*dly
	#last row and col are undefined. Replicate last value form previous row/col
	cellgr[rows,,]=cellgr[(rows-1),,]
	cellgr[,cols,]=cellgr[,(cols-1),]
	cellArea=sqrt(cellgr[,,1]^2+cellgr[,,2]^2+cellgr[,,3]^2)
	if (cArea) return(cellArea) else {
		for (i in 1:3) cellgr[,,i] = cellgr[,,i]/cellArea
		return(cellgr)  
	}
}
