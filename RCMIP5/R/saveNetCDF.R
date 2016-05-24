#' Save a cmip5data object to NetCDF format
#' 
#' There are at least three ways to save a \code{\link{cmip5data}} object.
#' First, \link{save} it. Second, use \link{as.data.frame}. Third, this function
#' will write out a new NetCDF file readable by any NetCDF-aware software.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param file Filename desired. If omitted one will be generated automatically.
#' @param path File path.
#' @param verbose logical. Print info as we go?
#' @param saveProvenance Save the provenance separately?
#' @param originalNames logical. Use original dimension names from file?
#' @return The fully-qualified filename that was written (invisible).
#' @details If no filename is provided, a meaningful one will be assigned based on the
#' CMIP5 naming convention (but appending 'RCMIP5'). \code{\link{loadCMIP5}} should be
#' able to read this file. If \code{saveProvenance} is specified, the provenance is saved
#' separately in a comma-separated file of the same name but appending "_prov.csv".
#' (Provenance messages are always saved as NetCDF file attributes.)
saveNetCDF <- function(x, file=NULL, path="./", verbose=FALSE, saveProvenance=TRUE, originalNames=FALSE) {
    
    # Sanity checks - class and length of parameters
    stopifnot(class(x)=="cmip5data")
    stopifnot(is.null(file) | (length(file)==1 & is.character(file)))
    stopifnot(length(path)==1 & is.character(path))
    stopifnot(length(verbose)==1 & is.logical(verbose))
    
    if(originalNames) 
        dimNames <- x$dimNames
    else
        dimNames <- c("lon", "lat", "Z", "time")
    if(verbose) cat("Writing with names", dimNames, "\n")
    
    # loadEnsemble() handles both ncdf and ncdf4, but saveNetCDF only supports the latter
    if(!require(ncdf4, quietly=!verbose)) {
        stop("This requires the 'ncdf4' package. Download it from CRAN or http://cirrus.ucsd.edu/~pierce/ncdf/")            
    }
    
    # Create meaningful filename, if necessary
    if(is.null(file)) {
        ensembles <- paste(x$ensembles, collapse="")
        pretty <- function(x) {formatC(round(x, 0), width=2, flag="0")}
        mintime <- paste0( floor(min(x$time)), pretty(min(x$time) %% 1 * 12 + 0.5))
        maxtime <- paste0( floor(max(x$time)), pretty(max(x$time) %% 1 * 12 + 0.5))
        file <- paste(x$variable, x$domain, x$model, x$experiment, ensembles, 
                      paste(mintime, maxtime, sep="-"), "RCMIP5.nc", sep="_")
    }
    fqfn <- paste(path, file, sep="/")
    
    # Define spatial dimensions, if present
    if(verbose) cat("Defining NetCDF dimensions...")
    dimlist <- list()
    if(!is.null(x$lon) & !is.null(x$lat)) {
        londim <- ncdf4::ncdim_def("lon", x$debug$lonUnit, x$lon)
        latdim <- ncdf4::ncdim_def("lat", x$debug$latUnit, x$lat)
        dimlist <- list(londim, latdim) # for now assume no Z/time        
    }
    
    # Define Z and time dimensions, if present
    if(!is.null(x$Z)) {
        Zdim <- ncdf4::ncdim_def(dimNames[3], x$debug$ZUnit, x$Z)
        dimlist <- list(londim, latdim, Zdim)
    }
    if(!is.null(x$time)) {
        timedim <- ncdf4::ncdim_def(dimNames[length(dimNames)], x$debug$timeUnit, x$debug$timeRaw, calendar=x$debug$calendarStr)
        dimlist[[length(dimlist)+1]] <- timedim     
    }
    
    if(verbose) cat(length(dimlist), "dimensions for", x$variable, "\n")
    
    # Define mandatory variable
    if(verbose) cat("Defining main NetCDF variable\n")
    valvar <- ncdf4::ncvar_def(x$variable, x$valUnit, dimlist)
    
    # Create the file and write mandatory variable
    # Note we make sure data is sorted correctly first
    if(verbose) cat("Creating and writing", file, "\n")
    nc <- ncdf4::nc_create(fqfn, valvar)    
    ncdf4::ncvar_put(nc, valvar, as.array(x))
    
    # Write spatial dimensions, if present
    if(!is.null(x$lon) & !is.null(x$lat)) {
        if(verbose) cat("Writing lon and lat\n")
        lonvar <- ncdf4::ncvar_def("lon", x$debug$lonUnit, londim)
        latvar <- ncdf4::ncvar_def("lat", x$debug$latUnit, londim)
        ncdf4::ncvar_put(nc, lonvar, x$lon)
        ncdf4::ncvar_put(nc, latvar, x$lat)        
    }
    
    # Write Z and time dimensions, if present
    if(!is.null(x$Z)) {
        if(verbose) cat("Writing Z\n")
        Zvar <- ncdf4::ncvar_def(dimNames[3], x$debug$ZUnit, Zdim)
        ncdf4::ncvar_put(nc, Zvar, x$Z) 
    }
    if(!is.null(x$time)) {
        if(verbose) cat("Writing time\n")
        timevar <- ncdf4::ncvar_def(dimNames[4], x$debug$timeUnit, timedim)
        ncdf4::ncvar_put(nc, timevar, x$time) 
    }
    
    # Get package version number, allowing that there might not be one
    pkgv <- "???"
    try({ pkgv <- packageVersion("RCMIP5") }, silent=T)
    
    # Write attributes
    if(verbose) cat("Writing attributes\n")    
    ncdf4::ncatt_put(nc, 0, "software", paste("Written by RCMIP5", pkgv, 
                                        "under", R.version.string, date()))
    if(!is.null(x$time)) {
        ncdf4::ncatt_put(nc, 0, "frequency", x$debug$timeFreqStr)
    }
    for(i in 1:nrow(x$provenance)) {
        ncdf4::ncatt_put(nc, 0, paste0("provenance", i), x$provenance[i, "message"])
    }
    
    # All done. Close file, update provenance, return filename
    ncdf4::nc_close(nc)
    if(verbose) cat("Wrote", round(file.info(fqfn)$size/1024/1024, 2), "MB\n")
    
    if(saveProvenance) {
        fqfn_prov <- paste0(fqfn, "_prov.csv")
        write.csv(as.data.frame(x$provenance), fqfn_prov)
        if(verbose) cat("Wrote provenance to", fqfn_prov, "\n")   
    }
    
    invisible(fqfn)
} # saveNetCDF
