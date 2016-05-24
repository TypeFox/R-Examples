doshade <- function(dem,sv,dl=0,sombra=dem) {
	if (!is.loaded('doshade')) {
		dyn.load('doshade.so')
	}
	if (nargs() < 2) {
		cat("USAGE: doshade(dem,sunvector,dl) \n")
		return()
	}
	switchdem=0
		if (class(dem)=="RasterLayer") {
			switchdem=1
			dproj=projection(dem)
			dext=extent(dem)
			dl = res(dem)[1]
			dem=raster::as.matrix(dem)
		}
	cols=ncol(dem)
	rows=nrow(dem)
	if (dl == 0){
			cat("Input data is not a RasterLayer, then I need the DEM resolution dl \n")
		return()
	}
	dem[is.na(dem)]=-999
	out = .Fortran("doshade",	
	dem=as.numeric(t(dem)),
	sunvector=as.vector(sv),
	cols=as.integer(cols),
	rows=as.integer(rows),
	dl=as.double(dl),
	sombra=as.numeric(dem),
	PACKAGE="insol")
	sombra=t(matrix(out$sombra,nrow=cols))
	if (switchdem){
		sombra=raster(sombra,crs=dproj)
		extent(sombra)=dext
	}
	sombra[sombra== -999]=NA
	return(sombra)
}
