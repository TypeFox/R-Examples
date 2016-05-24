#' Load a unique CMIP5 ensemble
#'
#' Loads the data for a particular CMIP5 experiment-variable-model-ensemble
#' combination (one or more files). Returns NULL and a warning if nothing matches.
#'
#' @param variable CMIP5 variable to load (required)
#' @param model CMIP5 model to load (required)
#' @param experiment CMIP5 experiment to load (required)
#' @param ensemble CMIP5 ensemble to load (required)
#' @param domain optinal CMIP5 domain to load (required)
#' @param path optional root of directory tree
#' @param recursive logical. Recurse into directories?
#' @param verbose logical. Print info as we go?
#' @param force.ncdf Force use of the less-desirable ncdf package for testing?
#' @param yearRange numeric of length 2. If supplied, load only these years of data inclusively between these years.
#' @return A \code{\link{cmip5data}} object, or \code{NULL} if nothing loaded
#' @details This function is the core of RCMIP5's data-loading. It loads all files matching
#' the experiment, variable, model, ensemble, and domain supplied by the caller.
#' @note This function is not intended to be called directly by the user; it will return
#' the \code{val} component as a multidimensional array, not a data frame.
#' @note This is an internal RCMIP5 function and not exported.
#' @keywords internal
loadEnsemble <- function(variable, model, experiment, ensemble, domain,
                         path='.', recursive=TRUE, verbose=FALSE, force.ncdf=FALSE,
                         yearRange=NULL) {
    
    # Sanity checks - make sure all parameters are correct class and length
    stopifnot(length(variable)==1 & is.character(variable))
    stopifnot(length(model)==1 & is.character(model))
    stopifnot(length(experiment)==1 & is.character(experiment))
    stopifnot(length(ensemble)==1 & is.character(ensemble))
    stopifnot(length(domain)==1 & is.character(domain))
    stopifnot(length(path)==1 & is.character(path)) # valid path?
    stopifnot(file.exists(path))
    stopifnot(length(recursive)==1 & is.logical(recursive))
    stopifnot(length(verbose)==1 & is.logical(verbose))
    stopifnot(length(force.ncdf)==1 & is.logical(force.ncdf))
    stopifnot(is.null(yearRange) | length(yearRange)==2 & is.numeric(yearRange))
    
    # We prefer to use the 'ncdf4' package, but Windows has problems with this,
    # ...so if it's not installed can also use 'ncdf'
    if(force.ncdf | !require(ncdf4, quietly=!verbose)) {
        if(require(ncdf, quietly=!verbose)) {
            # The ncdf and ncdf4 functions are mostly parameter-identical.
            # ...This makes things easy: we redefine the ncdf4 function
            # ...names to their ncdf equivalents
            .nc_open <- ncdf::open.ncdf
            .ncatt_get <- ncdf::att.get.ncdf
            .ncvar_get <- ncdf::get.var.ncdf
            .nc_close <- ncdf::close.ncdf
        } else {
            stop("No NetCDF (either 'ncdf4' or 'ncdf') package is available")
        }
    } else {
        .nc_open <- ncdf4::nc_open
        .ncatt_get <- ncdf4::ncatt_get
        .ncvar_get <- ncdf4::ncvar_get
        .nc_close <- ncdf4::nc_close
    }
    
    # List all files that match specifications
    fileList <- list.files(path=path, full.names=TRUE, recursive=recursive)
    
    # Match file names with valid CMIP5 patterns:
    # ...variable_domain_model_experiment_ensemble followed by either
    # ...a '_' or a '.' depending on whether a time period is specified or not.
    # ...Note that the '[_\\.]' match differentiates between
    # ...ensembles like 'r1i1p1' and 'r1i1p11', an unlikely but possible case.
    fileList <- fileList[grepl(pattern=sprintf('^%s_%s_%s_%s_%s[_\\.].*nc$',
                                               variable, domain, model,
                                               experiment, ensemble),
                               basename(fileList))]
    
    if(length(fileList)==0) {
        warning(paste("Could not find any matching files for",
                      variable, domain, model, experiment, ensemble))
        return(NULL)
    }
    
    # Get the domains of all files we want to load. Recall CMIP5 naming
    # ...conventions:variable_domain_model_experiment_ensemble_time.nc
    # ...
    domainCheck <- unname(vapply(unlist(fileList),
                                 function(x) { unlist(strsplit(basename(x), '_'))[2] },
                                 FUN.VALUE=''))
    
    # Check that we are only loading one domain. We check this before checking
    # other CMIP5 specifications because 'fx' domains will split on '_' to a
    # different number of strings then temporal domains.
    if(length(unique(domainCheck)) > 1) {
        stop('Domain is not unique: [', paste(unique(domainCheck), collapse=' '), ']\n')
    }
    
    # Get the number of pieces of CMIP5 information strings
    numSplits <- length(unlist(strsplit(basename(fileList[1]), '_')))
    
    # Split all file names based on '_' or '.'
    cmipName <- unname(vapply(unlist(fileList),
                              function(x) { unlist(strsplit(basename(x), '[_\\.]')) },
                              FUN.VALUE=rep('', length=numSplits+1)))
    
    # List what order the files appear in the name, for CMIP5 this will be:
    # ...variable_domain_model_experiment_ensemble_time.nc Note that we
    # ...aren't interested in checking the time string
    checkField <- list(variable=1, domain=2, model=3, experiment=4, ensemble=5)
    
    # Go through and make sure that we have unique variable, domain, model,
    # ...experiment, and ensemble strings for the set of files we are trying
    # ...to load. We want to avoid trying to load a file from ModelA and ModelB
    # ...as one cmip5data structure here. Also set the appropreate variables
    # ...so we can identify the cmip5data object correctly.
    for(checkStr in names(checkField)) {
        # Pull all unique strings
        tempStr <- unique(cmipName[checkField[[checkStr]],])
        if(length(tempStr) > 1) { # there should only be one!
            stop('[',checkStr, '] is not unique: [', paste(tempStr, collapse=' '), ']\n')
        } else {
            # Set a variable by the field name to the correct string for
            # ...this ensemble. For example there is now a variable
            # ...'variable' in the workspace containing the correct string
            # ...extracted from the first position of the file name.
            eval(parse(text=paste(checkStr, ' <- "', tempStr[1], '"', sep='')))
        }
    }
    
    # Go through and load the data
    val <- c() # variable to temporarily holds main data
    timeRaw <- c()
    timeArr <- c()
    ZUnit <- NULL
    valUnit <- NULL
    loadedFiles <- c()
    
    # Note that list.files returns a sorted list so these file should already
    # be in temporal order if the ensemble is split over multiple files.
    for(fileStr in fileList) {
        if(verbose) cat('Loading', fileStr, "\n")
        nc <- .nc_open(fileStr, write=FALSE)
        
        # Get variable names available
        varNames <- unlist(lapply(nc$var, FUN=function(x) { x$name }))
        
        # Get dimension names for 'variable'
        dimNames <- unlist(lapply(nc$var[[variable]]$dim, FUN=function(x) { x$name }))
        if(verbose) cat("-", variable, "dimension names:", dimNames, "\n")
        stopifnot(length(dimNames) %in% c(1, 2, 3, 4)) # that's all we know
        
        # Most, but not all, files have longitude and latitude. Load if available.
        lonArr <- NULL
        lonUnit <- NULL
        latArr <- NULL
        latUnit <- NULL
        if(length(dimNames) > 1) {
            # If 'lon' and 'lat' are available, use those. Otherwise load
            # whatever is specified by the main variable's dimensionality
            # TODO: this is a temporary (?) hack re issue #96
            if('lon' %in% varNames) dimNames[1] <- 'lon'
            if('lat' %in% varNames) dimNames[2] <- 'lat'
            
            lonArr <- .ncvar_get(nc, varid=dimNames[1])
            lonUnit <- .ncatt_get(nc, dimNames[1], 'units')$value
            latArr <- .ncvar_get(nc, varid=dimNames[2])
            latUnit <- .ncatt_get(nc, dimNames[2], 'units')$value
            
            # Some models provide two-dimensional arrays of their lon and lat values.
            # (Looking at you, GFDL.) If this occurs, strip down to 1
            if(length(dim(lonArr)) > 1) lonArr <- as.vector(lonArr[,1])
            if(length(dim(latArr)) > 1) latArr <- as.vector(latArr[1,])
        }
        
        # Get the time frequency. Note that this should be related to
        # ...the domain (ie 'mon' should be the frequency of the domain 'Amon').
        # ...In theory we could extract this from the domain
        timeFreqStr <- .ncatt_get(nc, varid=0, "frequency")$value
        
        # Non-fixed files have a time dimension to deal with:
        if(! timeFreqStr %in% 'fx') {
            # Get the time unit (e.g. 'days since 1860')
            timeName <- dimNames[length(dimNames)]
            timeUnit <- .ncatt_get(nc, timeName, 'units')$value
            # Get the type of calendar used (e.g. 'noleap')
            calendarStr <- .ncatt_get(nc, timeName, 'calendar')$value
            calendarUnitsStr <- .ncatt_get(nc, timeName, 'units')$value
            
            # Extract the number of days in a year
            if(grepl('^[^\\d]*\\d{3}[^\\d]day', calendarStr)) {
                calendarDayLength <- as.numeric(regmatches(calendarStr,
                                                           regexpr('\\d{3}', calendarStr)))
            } else {
                calendarDayLength <- 365
            }
            
            # Extract the year we are references in calendar
            # Set the default to year 1, month 1, day 1, hour 0, min 0, sec 0
            defaultCalendarArr <- c(1, 1, 1, 0, 0, 0)
            
            # Split the calandar unit string based on a '-', space, or ':'
            # ...this allows us to deal with YYYY-MM-DD hh:mm:ss, YYYY-MM-DD, or
            # ...YYYY-M-D
            calendarArr <- unlist(strsplit(calendarUnitsStr, split='[- :]'))
            
            # Check that the time is going to be in days otherwise latter
            # ...calculations for time array break
            stopifnot(any(grepl('day', calendarArr)))
            
            # extract just the digits
            calendarArr <- as.numeric(calendarArr[grepl('^\\d+$', calendarArr)])
            
            # swap the default values with the extracted values, we assume
            # ... that years are listed before months, before days, and so on
            temp <- defaultCalendarArr
            temp[1:length(calendarArr)] <- calendarArr
            calendarArr <- temp
            
            # calculate the decimal starting year
            startYr <- sum((calendarArr - c(0, 1, 1, 0, 0, 0))
                           / c(1, 12, calendarDayLength, calendarDayLength*24,
                               calendarDayLength*24*60, calendarDayLength*24*60*60))
            
            # Load the actual time
            thisTimeRaw <- .ncvar_get(nc, varid=timeName)
            # convert from days (we assume the units are days) to years
            thisTimeArr <- thisTimeRaw / calendarDayLength + startYr
        } else { # this is a fx variable. Set most things to NULL
            startYr <- NULL
            timeArr <- NULL
            timeUnit <- NULL
            calendarStr <- NULL
            calendarDayLength <- NULL
            calendarUnitsStr <- NULL
            dim(val) <- dim(val)[1:2]
            thisTimeRaw <- NULL
            thisTimeArr <- NULL
        }
        
        # Load the 4th dimension, if present:
        ZArr <- NULL
        if(length(dimNames) == 4) {
            ZArr <- .ncvar_get(nc, varid=dimNames[3])
            ZUnit <- .ncatt_get(nc, dimNames[3], 'units')$value
        }
        
        # If yearRange supplied, calculate filter for the data load below
        start <- NA
        count <- NA
        if(!is.null(yearRange) & !is.null(thisTimeArr)) {
            # User has requested to load a temporal subset of the data.
            # First question: does this file overlap at all?
            if(min(yearRange) > max(floor(thisTimeArr)) |
                   max(yearRange) < min(floor(thisTimeArr))) {
                if(verbose) cat("- skipping file because not in yearRange\n")
                next
            }
            
            # Calculate what positions in time array fall within yearRange
            # find first time match
            tstart <- match(min(yearRange), floor(thisTimeArr))
            if(is.na(tstart)) tstart <- 1
            # find last time match
            tend <- match(max(yearRange)+1, floor(thisTimeArr)) - 1
            if(is.na(tend)) tend <- length(thisTimeArr)
            
            # Construct the 'start' and 'count' arrays for ncvar_get below
            # (See ncvar_get documentation for what these mean.)
            ndims <- nc$var[[variable]]$ndims
            start <- c(rep(1, ndims-1), tstart)
            count <- c(rep(-1, ndims-1), tend-tstart+1)
            if(verbose) cat("- loading only timeslices", tstart, "-",tend, "\n")
            
            # Trim the already-loaded time arrays to match
            thisTimeArr <- thisTimeArr[tstart:tend]
            thisTimeRaw <- thisTimeRaw[tstart:tend]
        } # if year range
        
        # Update running time data
        if(!is.null(thisTimeRaw)) {
            timeRaw <- c(timeRaw, thisTimeRaw)
            timeArr <- c(timeArr, thisTimeRaw / calendarDayLength + startYr)
        }
        
        # Finally, load the actual data and its units
        vardata <- .ncvar_get(nc, varid=variable, start=start, count=count)
        if(verbose) cat("- data", dim(vardata), "\n")
        valUnit <- .ncatt_get(nc, variable, 'units')$value  # load units
        loadedFiles <- c(loadedFiles, basename(fileStr))
        
        # Restore any 'missing' dimensions (because not present, or length=1),
        # inserting NAs into dimNames to mark what wasn't present in file
        temp <- restoreMissingDims(dim(vardata), dimNames, lonArr, latArr, ZArr,
                                   thisTimeRaw, verbose)
        vardata <- array(vardata, dim=temp[["dims"]])
        dimNames <- temp[["dimNames"]]
        
        # Test that spatial dimensions are identical across files
        if(length(val) > 0 & length(dimNames) > 2) {
            stopifnot(all(dim(val)[1:(length(dim(val))-1)] ==
                              dim(vardata)[1:(length(dim(vardata))-1)]))
        }
        
        # Bind the main variable along time dimension to previously loaded data
        # Note that the time dimension, if present, is guaranteed to be last
        # ...see ncdf4 documentation
        val <- abind(val, vardata, along=length(dim(vardata)))
        
        .nc_close(nc)
    } # for filenames
    
    # If nothing loaded...
    if(length(val) == 0) return(NULL)
    
    x <- cmip5data(list(files=loadedFiles, val=unname(val), valUnit=valUnit,
                        lat=latArr, lon=lonArr, Z=ZArr, time=timeArr,
                        variable=variable, model=model, domain=domain,
                        experiment=experiment, ensembles=ensemble,
                        dimNames=dimNames,
                        
                        debug=list(startYr=startYr,
                                   lonUnit=lonUnit, latUnit=latUnit,
                                   ZUnit=ZUnit,
                                   timeUnit=timeUnit,
                                   timeFreqStr=timeFreqStr,
                                   calendarUnitsStr=calendarUnitsStr,
                                   calendarStr=calendarStr, timeRaw=timeRaw,
                                   calendarDayLength=calendarDayLength)
    ))
    
    # Add the provenance information, with a line for each loaded file
    for(f in fileList) {
        x <- addProvenance(x, paste("Loaded", basename(f)))
    }
    
    x
} # loadEnsemble


#' Restore missing and/or degenerate dimensions in the data
#'
#' @param dims the data array just loaded from the NetCDF
#' @param dimNames vector of dimensions names present in file
#' @param lonArr numeric vector of longitude values
#' @param latArr numeric vector of latitude values
#' @param ZArr numeric vector of Z values
#' @param thisTimeRaw numeric vector of time values
#' @param verbose logical. Print info as we go?
#' @return The data array with restored dimensions.
#' @note There are two cases to consider here. (1) If we load a dimension with only
#' one value (one month, one depth, etc) then that dimension will be dropped by
#' the .ncvar_get function (there's a 'collapse_degen' option available in ncdf4,
#' but not in ncdf). (2) A dimension is missing entirely, e.g. in a time-only file,
#' or a space-only grid area file. In either case, we want those dimensions (of
#' length 1) back in the data array.
#' @note This is an internal RCMIP5 function and not exported.
#' @keywords internal
restoreMissingDims <- function(dims, dimNames, lonArr, latArr, ZArr, thisTimeRaw, verbose) {
    if(is.null(lonArr) | length(lonArr) == 1 ) {
        if(verbose) cat("- restoring dimension for lon\n")
        dims <- c(1, dims)
    }
    if(is.null(latArr) | length(latArr) == 1 ) {
        if(verbose) cat("- restoring dimension for lat\n")
        dims <- c(dims[1], 1, dims[2:length(dims)])
    }
    if(length(ZArr) == 1) {
        if(verbose) cat("- restoring dimension for Z\n")
        if(length(dims) >= 3)
            dims <- c(dims[1:2], 1, dims[3:length(dims)])
        else
            dims <- c(dims, 1)
    }
    if(length(thisTimeRaw) == 1) {
        if(verbose) cat("- adding extra dimension for time\n")
        dims <- c(dims, 1)
    }
    
    # At this point, we've restored all dimensions dropped due to length 1 issues
    # But we want all data moving through RCMIP5 to have four dimensions
    if(length(dims) == 1) {  # assume time only
        if(verbose) cat("- adding extra dimensions for lon, lat, Z\n")
        dims <- c(1, 1, 1, dims)
        dimNames <- c(NA, NA, NA, dimNames)
    } else if(length(dims) == 2) { # assume lon, lat
        if(verbose) cat("- adding extra dimensions for Z, time\n")
        dims <- c(dims, 1, 1)
        dimNames <- c(dimNames, NA, NA)
    } else if(length(dims) == 3) { # assume lon, lat, time
        if(verbose) cat("- adding extra dimension for Z\n")
        dims <- c(dims[1:2], 1, dims[3])
        dimNames <- c(dimNames[1:2], NA, dimNames[3])
    } else if(length(dims) == 4) { # assume lon, lat, Z, time
        # no change needed
    } else
        stop("Variable dimensions out of bounds!")
    
    list(dims=dims, dimNames=dimNames)
} # restoreMissingDimensions
