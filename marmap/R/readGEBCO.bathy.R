readGEBCO.bathy <- function(file, resolution=1){

	# require(ncdf)

	# check resolution value ## is.wholenumber function from is.integer {base}
	"Argument 'resolution' must be a positive integer\n" -> message
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol	
	# if resolution is <1 OR not a whole number:
	if(resolution<1 | !is.wholenumber(resolution)) stop(message) else {

		# get data from netCDF file
		nc <- ncdf4::nc_open(file)
		# ncells <- length(ncdf::get.var.ncdf(nc, "xysize"))
		z <- ncdf4::ncvar_get(nc,"elevation")
		lat <- ncdf4::ncvar_get(nc, "lat")
		lon <- ncdf4::ncvar_get(nc, "lon")
	
		# dimensions of the matrix, depending on type of database
		# if(db == "GEBCO_1min") db.scale <- 1/60
		# if(db == "GEBCO_08")   db.scale <- 1/120
		# xdim <- seq(xrg[1],xrg[2], by=db.scale)
		# ydim <- seq(yrg[1],yrg[2], by=db.scale)
		
		# # for some reason sometimes the z vector is shorter than the product of the
		# # length of the latitude and longitude vectors
		# # so we check that z and xy are compatible, otherwise we crop the matrix dimentions.
		# if(length(xdim) * length(ydim) != ncells) {
		# 	warning("The number of cells in the .nc file (", ncells , ") do not correspond exactly to the range of latitude x longitude values (", length(xdim) * length(ydim), ")...\n  Cropping the matrix dimentions (see ?readGEBCO.bathy for details)...")
		# 	xdim<-xdim[-length(xdim)];ydim<-ydim[-length(ydim)] # removing last x and y values
		# }

		# build matrix 
		mat <- matrix(data=rev(z), nrow=length(lon),ncol=length(lat), byrow=F, dimnames=list(rev(lon),rev(lat)))
		# mat <- mat[,ncol(mat):1] # (merci benoit!)
		# rownames(mat) <- xdim
		# colnames(mat) <- ydim
		mat <- check.bathy(mat[,ncol(mat):1])
	
		# adapt the resolution
		xindex <- seq(1,length(lon),by=resolution)
		yindex <- seq(1,length(lat),by=resolution)
		mat[xindex,yindex] -> final
		
		class(final) <- "bathy"
		return(final)
	}
}