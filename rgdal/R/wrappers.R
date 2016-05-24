

# thin wrappers 121001 Robert J. Hijmans

.gd_SetNoDataValue <- function(object, NAflag) {
	.Call("RGDAL_SetNoDataValue", object, as.double(NAflag), PACKAGE="rgdal")
}


.gd_SetGeoTransform <- function(object, geotrans) {
   .Call("RGDAL_SetGeoTransform", object, geotrans, PACKAGE="rgdal")
}


.gd_SetProject <- function(object, proj4string) {
    .Call("RGDAL_SetProject", object, proj4string, PACKAGE="rgdal")
}


.gd_SetStatistics <- function(object, statistics) {
	.Call("RGDAL_SetStatistics", object, as.double(statistics), PACKAGE="rgdal")
}

.gd_SetRasterColorTable <- function(object, icT) {
    if (!is.matrix(icT)) {
        stopifnot(is.character(icT))
        icT <- t(col2rgb(icT))
    }
    stopifnot(is.matrix(icT))
    stopifnot(storage.mode(icT) == "integer")
    ricT <- nrow(icT)
    cicT <- ncol(icT)
    stopifnot(cicT == 3 || cicT == 4)
    .Call("RGDAL_SetRasterColorTable", object, icT, ricT, cicT,
        PACKAGE = "rgdal")
}

.gd_SetCategoryNames <- function(object, icN) {
    stopifnot(is.character(icN))
    .Call("RGDAL_SetCategoryNames", object, icN, PACKAGE = "rgdal")
}


GDALcall <- function(object, option, ...) {
	if (option == 'SetNoDataValue') {
		.gd_SetNoDataValue(object, ...)
	} else if (option == 'SetGeoTransform') {
		.gd_SetGeoTransform(object, ...)
	} else if (option == 'SetProject') {
		.gd_SetProject(object, ...)
	} else if (option == 'SetStatistics') {
		.gd_SetStatistics(object, ...)
	} else if (option == 'SetRasterColorTable') {
		.gd_SetRasterColorTable(object, ...)
	} else if (option == 'SetCategoryNames') {
		.gd_SetCategoryNames(object, ...)
	} else {
		stop('invalid option')
	}
}




.gd_transform <- function(projfrom, projto, n, x, y, z=NULL) {
  if (!get("has_proj_def.dat", envir=.RGDAL_CACHE)) {
      projfrom <- proj_def_bug_fix(projfrom)
      projto <- proj_def_bug_fix(projto)
  }
  if (is.null(z)) .Call("transform", projfrom, projto, n, x, y, NULL, PACKAGE="rgdal")
  else .Call("transform", projfrom, projto, n, x, y, z, PACKAGE="rgdal")
}

# exported version
rawTransform <- function(projfrom, projto, n, x, y, z=NULL) {
        if (!get("has_proj_def.dat", envir=.RGDAL_CACHE)) {
            projfrom <- proj_def_bug_fix(projfrom)
            projto <- proj_def_bug_fix(projto)
        }
	if (is.null(z)) .Call("transform", projfrom, projto, n, x, y, NULL, PACKAGE="rgdal")
	else .Call("transform", projfrom, projto, n, x, y, z, PACKAGE="rgdal")
}

