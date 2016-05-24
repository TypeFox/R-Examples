#===============================================================================#
#										#
#  Name:       RNetCDF.R							#
#										#
#  Version:    1.8-2								#
#										#
#  Purpose:    NetCDF interface for R.						#
#										#
#  Author:     Pavel Michna (michna@giub.unibe.ch)				#
#              Milton Woods (m.woods@bom.gov.au)                                #
#										#
#  Copyright:  (C) 2004-2014 Pavel Michna					#
#										#
#===============================================================================#
#										#
#  This program is free software; you can redistribute it and/or modify 	#
#  it under the terms of the GNU General Public License as published by 	#
#  the Free Software Foundation; either version 2 of the License, or		#
#  (at your option) any later version.						#
#										#
#  This program is distributed in the hope that it will be useful,		#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of		#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		#
#  GNU General Public License for more details. 				#
#										#
#  You should have received a copy of the GNU General Public License		#
#  along with this program; if not, write to the Free Software			#
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA	#
#										#
#===============================================================================#
#  Implementation and Revisions 						#
#-------------------------------------------------------------------------------#
#  Author   Date       Description						#
#  ------   ----       -----------						#
#  pm       12/06/04   First implementation					#
#  pm       09/07/04   Support scalar variables	    		                #
#  pm       21/07/04   Changed error handling					#
#  pm       28/07/04   Minor modifications					#
#  pm       12/09/04   New na.mode=3 and collapse=TRUE/FALSE in var.get.nc()	#
#  pm       24/07/06   Handling dates in string form (udunits)           	#
#  mw       14/04/08   Added new modes (large, prefill, share)                  #
#                      to nc_open and nc_create                                 #
#  pm       24/11/10   Added new option enddef to att and dim/var definitions   #
#  pm       01/12/10   Removed option enddef, checking in C code for mode       #
#  pm       14/02/12   Corrected bug in att.delete.nc                           #
#  pm       02/06/12   Added function read.nc()                                 #
#  pm       16/07/12   Added packing/unpacking of data (code from mw)           #
#  mw       21/08/14   Allow reading of character vector or scalar              #
#  mw       05/09/14   Support reading and writing raw character arrays         #
#  mw       08/09/14   Handle reading and writing of zero-sized arrays          #
#  mw       24/01/16   Support conversion of timestamps to/from POSIXct         #
#										#
#===============================================================================#


#===============================================================================#
#  NetCDF library functions							#
#===============================================================================#

#-------------------------------------------------------------------------------#
#  att.copy.nc()                                                                #
#-------------------------------------------------------------------------------#

att.copy.nc <- function(ncfile.in, variable.in, attribute, ncfile.out,
                        variable.out)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile.in ) == "NetCDF")
    stopifnot(class(ncfile.out) == "NetCDF")
    stopifnot(is.character(attribute)    || is.numeric(attribute)   )
    stopifnot(is.character(variable.in)  || is.numeric(variable.in) )
    stopifnot(is.character(variable.out) || is.numeric(variable.out))
    
    #-- Get the varids as integer if necessary, handle global attributes -------#
    varid.in     <- NULL
    varid.out    <- NULL
    globflag.in  <- 0
    globflag.out <- 0
    
    if(is.character(variable.in) && variable.in != "NC_GLOBAL")
        varid.in <- try(var.inq.nc(ncfile.in, variable.in)$id)
    else
        varid.in <- variable.in

    if(is.character(variable.in) && variable.in == "NC_GLOBAL") {
        globflag.in <-  1
        varid.in    <- -1
    }
    
    if(is.character(variable.out) && variable.out != "NC_GLOBAL")
        varid.out <- try(var.inq.nc(ncfile.out, variable.out)$id)
    else
        varid.out <- variable.out

    if(is.character(variable.out) && variable.out == "NC_GLOBAL") {
        globflag.out <-  1
        varid.out    <- -1
    }
	
    if(class(varid.in) == "try-error" || class(varid.out) == "try-error")
        return(invisible(NULL))
    if(is.null(varid.in) || is.null(varid.out))
        return(invisible(NULL))

    #-- Get the attribute name if necessary ------------------------------------#
    if(is.character(attribute))
        attname <- attribute
    else
        attname <- try(att.inq.nc(ncfile.in, variable.in, attribute)$name)
    
    if(class(attname) == "try-error" || is.null(attname))
        return(invisible(NULL))
    
    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_copy_att",
                as.integer(ncfile.in),
		as.integer(varid.in),
		as.integer(globflag.in),
		as.character(attname),
                as.integer(ncfile.out),
		as.integer(varid.out),
		as.integer(globflag.out),
		PACKAGE="RNetCDF")
		
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  att.delete.nc()                                                              #
#-------------------------------------------------------------------------------#

att.delete.nc <- function(ncfile, variable, attribute)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(variable)  || is.numeric(variable) )
    stopifnot(is.character(attribute) || is.numeric(attribute))
    
    #-- Get the varid as integer if necessary, handle global attributes --------#
    varid     <- NULL
    globflag  <- 0
    
    if(is.character(variable) && variable != "NC_GLOBAL")
        varid <- try(var.inq.nc(ncfile, variable)$id)
    else
        varid <- variable

    if(is.character(variable) && variable == "NC_GLOBAL") {
        globflag <-  1
        varid    <- -1
    }

    if(class(varid) == "try-error" || is.null(varid))
        return(invisible(NULL))

    #-- Get the attribute name if necessary ------------------------------------#
    if(is.character(attribute))
        attname <- attribute
    else
        attname <- try(att.inq.nc(ncfile, variable, attribute)$name)
    
    if(class(attname) == "try-error" || is.null(attname))
        return(invisible(NULL))

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_delete_att",
                as.integer(ncfile),
		as.integer(varid),
		as.integer(globflag),
		as.character(attname),
		PACKAGE="RNetCDF")
    
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  att.get.nc()                                                                 #
#-------------------------------------------------------------------------------#

att.get.nc <- function(ncfile, variable, attribute)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(variable)  || is.numeric(variable) )
    stopifnot(is.character(attribute) || is.numeric(attribute))
    
    #-- Get the varid as integer if necessary, handle global attributes --------#
    varid     <- NULL
    globflag  <- 0
    
    if(is.character(variable) && variable != "NC_GLOBAL")
        varid <- try(var.inq.nc(ncfile, variable)$id)
    else
        varid <- variable

    if(is.character(variable) && variable == "NC_GLOBAL") {
        globflag <-  1
        varid    <- -1
    }

    if(class(varid) == "try-error" || is.null(varid))
        return(invisible(NULL))

    #-- Inquire the attribute to get its name and storage mode -----------------#
    attinfo <- try(att.inq.nc(ncfile, variable, attribute))
    if(class(attinfo) == "try-error" || is.null(attinfo))
        return(invisible(NULL))

    ifelse(is.character(attribute),
        attname <- attribute, attname <- attinfo$name)
    
    ifelse(attinfo$type == "NC_CHAR", numflag <- 0, numflag <- 1);
    
    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_get_att",
	        as.integer(ncfile),
		as.integer(varid),
		as.character(attname),
		as.integer(numflag),
		as.integer(globflag),
		PACKAGE="RNetCDF")
    
    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        nc$status <- NULL
	nc$errmsg <- NULL
        return(nc$value)
    } else
	stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  att.inq.nc()                                                                 #
#-------------------------------------------------------------------------------#

att.inq.nc <- function(ncfile, variable, attribute)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(variable)  || is.numeric(variable) )
    stopifnot(is.character(attribute) || is.numeric(attribute))
    
    #-- Look if handle attribute by name or ID ---------------------------------#
    attid   <- -1
    attname <- ""
    ifelse(is.character(attribute), nameflag <- 1, nameflag <- 0)
    ifelse(is.character(attribute), attname <- attribute, attid <- attribute)

    #-- Get the varid as integer if necessary, handle global attributes --------#
    varid     <- NULL
    globflag  <- 0
    
    if(is.character(variable) && variable != "NC_GLOBAL")
        varid <- try(var.inq.nc(ncfile, variable)$id)
    else
        varid <- variable

    if(is.character(variable) && variable == "NC_GLOBAL") {
        globflag <-  1
        varid    <- -1
    }

    if(class(varid) == "try-error" || is.null(varid))
        return(invisible(NULL))

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_inq_att",
        	as.integer(ncfile),
		as.integer(varid),
		as.character(attname),
		as.integer(attid),
		as.integer(nameflag),
		as.integer(globflag),
		PACKAGE="RNetCDF")

    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        nc$status <- NULL
	nc$errmsg <- NULL
        return(nc)
    } else
	stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  att.put.nc()                                                                 #
#-------------------------------------------------------------------------------#

att.put.nc <- function(ncfile, variable, name, type, value)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(variable) || is.numeric(variable))
    stopifnot(is.character(name))
    stopifnot(is.character(type))
    stopifnot(is.character(value) || is.numeric(value))

    #-- Get the varids as integer if necessary, handle global attributes -------#
    varid     <- NULL
    globflag  <- 0
    
    if(is.character(variable) && variable != "NC_GLOBAL")
        varid <- try(var.inq.nc(ncfile, variable)$id)
    else
        varid <- variable

    if(is.character(variable) && variable == "NC_GLOBAL") {
        globflag <-  1
        varid    <- -1
    }

    if(class(varid) == "try-error" || is.null(varid))
        return(invisible(NULL))

    #-- Determine if attribute is numeric or character -------------------------#
    ifelse(is.numeric(value), numflag <- 1, numflag <- 0)

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_put_att",
	        as.integer(ncfile),
		as.integer(varid),
		as.character(name),
		as.character(type),
		as.integer(length(value)),
		as.integer(numflag),
		as.integer(globflag),
		value,
		PACKAGE="RNetCDF")
    
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  att.rename.nc()                                                              #
#-------------------------------------------------------------------------------#

att.rename.nc <- function(ncfile, variable, attribute, newname)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(variable)  || is.numeric(variable) )
    stopifnot(is.character(attribute) || is.numeric(attribute))
    stopifnot(is.character(newname))

    #-- Get the varid as integer if necessary, handle global attributes --------#
    varid     <- NULL
    globflag  <- 0
    
    if(is.character(variable) && variable != "NC_GLOBAL")
        varid <- try(var.inq.nc(ncfile, variable)$id)
    else
        varid <- variable

    if(is.character(variable) && variable == "NC_GLOBAL") {
        globflag <-  1
        varid    <- -1
    }
	
    if(class(varid) == "try-error" || is.null(varid))
        return(invisible(NULL))

    #-- Get the attribute name if necessary ------------------------------------#
    if(is.character(attribute))
        attname <- attribute
    else
        attname <- try(att.inq.nc(ncfile, variable, attribute)$name)
    
    if(class(attname) == "try-error" || is.null(attname))
        return(invisible(NULL))

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_rename_att",
                as.integer(ncfile),
		as.integer(varid),
		as.integer(globflag),
		as.character(attname),
                as.character(newname),
		PACKAGE="RNetCDF")
    
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  close.nc()                                                                   #
#-------------------------------------------------------------------------------#

close.nc <- function(con, ...)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(con) == "NetCDF") 

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_close",
                as.integer(con),
		PACKAGE="RNetCDF")

    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  create.nc()                                                                  #
#-------------------------------------------------------------------------------#

create.nc <- function(filename, clobber=TRUE, large=FALSE, share=FALSE,
                      prefill=TRUE)
{
    #-- Convert logical values to integers -------------------------------------#
    iclobber <- ifelse(clobber == TRUE, 1, 0)    ## Overwrite existing file (y/n)
    ilarge   <- ifelse(large   == TRUE, 1, 0)
    ishare   <- ifelse(share   == TRUE, 1, 0)
    iprefill <- ifelse(prefill == TRUE, 1, 0)

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_create",
		as.character(filename),
		as.integer(iclobber),
                as.integer(ilarge),
                as.integer(ishare),
                as.integer(iprefill),
		PACKAGE="RNetCDF")
	     
    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        ncfile <- nc$ncid
        attr(ncfile, "class") <- "NetCDF"
	return(ncfile)
    } else
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  dim.def.nc()                                                                 #
#-------------------------------------------------------------------------------#

dim.def.nc <- function(ncfile, dimname, dimlength=1, unlim=FALSE)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(dimname))
    stopifnot(is.numeric(dimlength))
    stopifnot(is.logical(unlim))

    #-- Handle UNLIMITED -------------------------------------------------------#
    unlimflag   <- ifelse(unlim  == TRUE, 1, 0)
    ncdimlength <- ifelse(unlim  == TRUE, 0, dimlength)
    
    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_def_dim",
		as.integer(ncfile),
		as.character(dimname),
		as.integer(ncdimlength),
		as.integer(unlimflag),
		PACKAGE="RNetCDF")

    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  dim.inq.nc()                                                                 #
#-------------------------------------------------------------------------------#

dim.inq.nc <- function(ncfile, dimension)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(dimension) || is.numeric(dimension))
    
    #-- Look if handle dimension by name or ID ---------------------------------#
    dimid   <- -1
    dimname <- ""
       
    ifelse(is.character(dimension), nameflag <- 1, nameflag <- 0)
    ifelse(is.character(dimension), dimname <- dimension, dimid <- dimension)
    
    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_inq_dim",
        	as.integer(ncfile),
		as.integer(dimid),
		as.character(dimname),
		as.integer(nameflag),
		PACKAGE="RNetCDF")
    
    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        nc$status <- NULL
	nc$errmsg <- NULL
	
	nc$unlim <- ifelse(nc$unlim == 1, TRUE, FALSE)
	
        return(nc)
    } else
	stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  dim.rename.nc()                                                              #
#-------------------------------------------------------------------------------#

dim.rename.nc <- function(ncfile, dimension, newname)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(dimension) || is.numeric(dimension))
    stopifnot(is.character(newname))

    #-- Look if handle dimension by name or ID ---------------------------------#
    dimid   <- -1
    dimname <- ""
       
    ifelse(is.character(dimension), nameflag <- 1, nameflag <- 0)
    ifelse(is.character(dimension), dimname <- dimension, dimid <- dimension)
    
    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_rename_dim",
        	as.integer(ncfile),
		as.integer(dimid),
		as.character(dimname),
		as.integer(nameflag),
		as.character(newname),
		PACKAGE="RNetCDF")

    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  file.inq.nc()                                                                #
#-------------------------------------------------------------------------------#

file.inq.nc <- function(ncfile)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_inq_file",
                as.integer(ncfile),
		PACKAGE="RNetCDF")

    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        nc$status <- NULL
	nc$errmsg <- NULL
        
	if(nc$unlimdimid == -1)
            nc$unlimdimid <- NA
	
	return(nc)
    } else
	stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  open.nc()                                                                    #
#-------------------------------------------------------------------------------#

open.nc <- function(con, write=FALSE, share=FALSE, prefill=TRUE, ...)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(is.character(con))
    stopifnot(is.logical(write))
    stopifnot(is.logical(share))
    stopifnot(is.logical(prefill))

    #-- Open read only (y/n) ---------------------------------------------------#
    iwrite   <- ifelse(write   == TRUE, 1, 0)
    ishare   <- ifelse(share   == TRUE, 1, 0)
    iprefill <- ifelse(prefill == TRUE, 1, 0)

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_open",
		as.character(con),
		as.integer(iwrite),
		as.integer(ishare),
		as.integer(iprefill),
		PACKAGE="RNetCDF")
	     
    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        ncfile <- nc$ncid
        attr(ncfile, "class") <- "NetCDF"
	return(ncfile)
    } else
       stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  print.nc()                                                                   #
#-------------------------------------------------------------------------------#

print.nc <- function(x, ...)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(x) == "NetCDF") 

    #-- Inquire about the file -------------------------------------------------#
    fileinfo <- try(file.inq.nc(x))
    if(class(fileinfo) == "try-error" || is.null(fileinfo))
        return(invisible(NULL))

    #-- Inquire about all dimensions -------------------------------------------#
    if(fileinfo$ndims != 0) {
	cat("dimensions:\n")
	for(i in 0:(fileinfo$ndims-1)) {
            diminfo <- dim.inq.nc(x, i)
            if(diminfo$unlim == FALSE)
		cat("        ", diminfo$name, " = ", diminfo$length, " ;\n",
		    sep="")
	    else
	        cat("        ", diminfo$name, " = UNLIMITED ; // (", 
		    diminfo$length, " currently)\n", sep="")
	}
    }
    
    #-- Inquire about all variables --------------------------------------------#
    if(fileinfo$nvars != 0) {
	cat("variables:\n")
	for(i in 0:(fileinfo$nvars-1)) {
            varinfo <- var.inq.nc(x, i)
	    vartype <- substring(tolower(varinfo$type), 4)
	    cat("        ", vartype, " ", varinfo$name, sep="")
	    if(varinfo$ndims > 0) {
	        cat("(")
		for(j in 1:(varinfo$ndims-1))
	            if(j > 0 && varinfo$ndims > 1)
			cat(dim.inq.nc(x, varinfo$dimids[j])$name, ", ", sep="")
		cat(dim.inq.nc(x, varinfo$dimids[varinfo$ndims])$name, sep="")
	        cat(")")
	    }
	    cat(" ;\n")
	    
	    #-- Inquire about variable attributes ------------------------------#
	    if(varinfo$natts != 0) {
	        for(j in 0:(varinfo$natts-1)) {
		    attinfo <- att.inq.nc(x, i, j)
		    cat(rep(" ", 16), varinfo$name, ":", attinfo$name, sep="")
		    if(attinfo$type == "NC_CHAR")
		        cat(" = \"", att.get.nc(x, i, j), "\" ;\n", sep="")
		    else
		        cat( " = ", att.get.nc(x, i, j), " ;\n", sep="")
		}
	    }
	}
    }

    #-- Inquire about gobal attributes------------------------------------------#
    if(fileinfo$ngatts != 0) {
        cat("\n// global attributes:\n")
        i <- "NC_GLOBAL"
	for(j in 0:(fileinfo$ngatts-1)) {
	    attinfo <- att.inq.nc(x, i, j)
	    cat(rep(" ", 16),  ":", attinfo$name, sep="")
	    if(attinfo$type == "NC_CHAR")
		cat(" = \"", att.get.nc(x, i, j), "\" ;\n", sep="")
	    else
		cat( " = ", att.get.nc(x, i, j), " ;\n", sep="")
	}
    }
}


#-------------------------------------------------------------------------------#
#  sync.nc()                                                                    #
#-------------------------------------------------------------------------------#

sync.nc <- function(ncfile)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF")

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_sync",
                as.integer(ncfile),
		PACKAGE="RNetCDF")
		
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  var.def.nc()                                                                 #
#-------------------------------------------------------------------------------#

var.def.nc <- function(ncfile, varname, vartype, dimensions)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF")
    stopifnot(is.character(varname))
    stopifnot(is.character(vartype))
    
    if(any(is.na(dimensions)) && length(dimensions) != 1)
        stop("NAs not allowed in dimensions unless defining a scalar variable",
	    call.=FALSE)

    if(!any(is.na(dimensions)))
	stopifnot(mode(dimensions) == "character" || 
	    mode(dimensions) == "numeric")
    
    #-- Determine dimids from dimname if necessary, handle scalar variables ----#
    dimids <- vector()
    if(any(is.na(dimensions))) {
        dimids <- -1
	ndims  <-  0
    } else {	
	if(mode(dimensions) == "numeric")
            dimids <- dimensions
	else
            for(i in seq_along(dimensions))
        	try(dimids[i] <- dim.inq.nc(ncfile, dimensions[i])$id,
	            silent=TRUE)

	if(length(dimids) != length(dimensions))
            stop("Could not determine all dimension ids", call.=FALSE)
	    
	dimids <- dimids[length(dimids):1]                   ## R to C convention
	ndims  <- length(dimids)
    }
    
    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_def_var",
        	as.integer(ncfile),
		as.character(varname),
		as.character(vartype),
		as.integer(ndims),
		as.integer(dimids),
		PACKAGE="RNetCDF")
		
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  var.get.nc()                                                                 #
#-------------------------------------------------------------------------------#

var.get.nc <- function(ncfile, variable, start=NA, count=NA, na.mode=0,
                       collapse=TRUE, unpack=FALSE, rawchar=FALSE)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF")
    stopifnot(is.character(variable) || is.numeric(variable))
    stopifnot(is.numeric(start) || is.logical(start))
    stopifnot(is.numeric(count) || is.logical(count))
    stopifnot(is.logical(collapse))
    stopifnot(is.logical(unpack))
    stopifnot(is.logical(rawchar))
    
    stopifnot(isTRUE(na.mode %in% c(0,1,2,3))) 

    #-- Inquire the variable ---------------------------------------------------#
    varinfo <- try(var.inq.nc(ncfile, variable))

    if(class(varinfo) == "try-error" || is.null(varinfo))
        return(invisible(NULL))
    
    ndims      <- varinfo$ndims

    #-- Get the varid as integer if necessary ----------------------------------#
    ifelse(is.character(variable), varid <- varinfo$id, varid <- variable)

    #-- Replace NA by defined limits -------------------------------------------#
    if(isTRUE(is.na(start))) {
      start <- rep(NA, ndims)
    }
    if(isTRUE(is.na(count))) {
      count <- rep(NA, ndims)
    }

    if(length(start) != ndims || length(count) != ndims)
      stop("Length of start/count is not ndims", call.=FALSE)

    start[is.na(start)] <- 1
    for (idim in seq_len(ndims)) {
      if (is.na(count[idim])) {
        count[idim] <- dim.inq.nc(ncfile, varinfo$dimids[idim])$length
      }
    }

    #-- Switch from R to C convention ------------------------------------------#
    c.start <- rev(start-1)
    c.count <- rev(count)

    #-- C function calls -------------------------------------------------------#
    if(varinfo$type == "NC_CHAR") {
        nc <- .Call("R_nc_get_vara_text",
        	    as.integer(ncfile),
		    as.integer(varid),
		    as.integer(c.start),
		    as.integer(c.count),
		    as.integer(ndims),
                    as.integer(rawchar),
		    PACKAGE="RNetCDF")
    } else {
	nc <- .Call("R_nc_get_vara_double",
        	    as.integer(ncfile),
		    as.integer(varid),
		    as.integer(c.start),
		    as.integer(c.count),
		    as.integer(ndims),
		    PACKAGE="RNetCDF")
    }

    #-- Adjust data ------------------------------------------------------------#
    if(nc$status == 0) {
	
	#-- Convert missing value to NA if defined in NetCDF file --------------#
        if (is.numeric(nc$data) && na.mode < 3) { 
	  tolerance <- 1*10^-5                              ## Allow rounding error
	  na.flag   <- 0

	  missval.flag <- 0
	  fillval.flag <- 0
	  
	  fillval <- try(att.inq.nc(ncfile, varinfo$name, "_FillValue"   ), 
	      silent=TRUE)
	  missval <- try(att.inq.nc(ncfile, varinfo$name, "missing_value"), 
	      silent=TRUE)

	  if(!(class(fillval) == "try-error"))
	      if(!is.null(fillval))
		  fillval.flag <- 1
	  if(!(class(missval) == "try-error"))
	      if(!is.null(missval))
		  missval.flag <- 1

	  if(na.mode == 0 && missval.flag == 1) {
	      na.value <- att.get.nc(ncfile, varinfo$name, "missing_value")
	      na.flag  <- 1
	  }
	  if(na.mode == 0 && fillval.flag == 1) {
	      na.value <- att.get.nc(ncfile, varinfo$name, "_FillValue")
	      na.flag  <- 1
	  }

	  if(na.mode == 1 && fillval.flag == 1) {
	      na.value <- att.get.nc(ncfile, varinfo$name, "_FillValue")
	      na.flag  <- 1
	  }

	  if(na.mode == 2 && missval.flag == 1) {
	      na.value <- att.get.nc(ncfile, varinfo$name, "missing_value")
	      na.flag  <- 1
	  }

	  if(na.flag == 1) {
	      nc$data[abs(nc$data-as.numeric(na.value)) < tolerance] <- NA 
	  }
        }

	#-- Unpack variables if requested (missing values are preserved) -------#
        if(unpack && is.numeric(nc$data)) {
            offset <- try(att.inq.nc(ncfile, varinfo$name, "add_offset"),
                          silent=TRUE)
            scale  <- try(att.inq.nc(ncfile, varinfo$name, "scale_factor"),
                          silent=TRUE)
            if((!inherits(offset, "try-error")) &&
               (!inherits(scale,  "try-error"))) {
	        add_offset   <- att.get.nc(ncfile, varinfo$name, "add_offset")
		scale_factor <- att.get.nc(ncfile, varinfo$name, "scale_factor")
                nc$data      <- nc$data*scale_factor+add_offset
            }
        }

	#-- Set dimensions, collapse degenerate dimensions ---------------------#
        if (is.character(nc$data) && ndims > 0) {
          # Drop string length dimension
          datadim <- count[-1]
        } else {
          datadim <- count
        }
        if (collapse) {
          # Drop singleton dimensions
          datadim <- datadim[datadim!=1]
        }
        if (length(datadim)<1) {
          # For compatibility with code written for RNetCDF<=1.6.x,
          # scalars and vectors always have a dimension attribute:
          datadim <- length(nc$data)
        }
        dim(nc$data) <- datadim

    #-- Return object if no error ----------------------------------------------#
	return(nc$data)
    } else
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  var.inq.nc()                                                                 #
#-------------------------------------------------------------------------------#

var.inq.nc <- function(ncfile, variable)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF")
    stopifnot(is.character(variable) || is.numeric(variable))

    #-- Look if handle variable by name or ID ----------------------------------#
    varid   <- -1
    varname <- ""

    ifelse(is.character(variable), nameflag <- 1, nameflag <- 0)
    ifelse(is.character(variable), varname <- variable, varid <- variable)

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_inq_var",
        	as.integer(ncfile),
		as.integer(varid),
		as.character(varname),
		as.integer(nameflag),
		PACKAGE="RNetCDF")

    #-- Return object if no error ----------------------------------------------#
    if(nc$status == 0) {
        nc$status <- NULL
	nc$errmsg <- NULL
        
	if(nc$ndims > 0)
            nc$dimids <- nc$dimids[(nc$ndims):1]             ## C to R convention
        else
            nc$dimids <- NA
	    
	return(nc)
    } else
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  var.put.nc()                                                                 #
#-------------------------------------------------------------------------------#

var.put.nc <- function(ncfile, variable, data, start=NA, count=NA, na.mode=0,
                       pack=FALSE)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF")
    stopifnot(is.character(variable) || is.numeric(variable))
    stopifnot(is.numeric(data) || is.character(data) || is.raw(data) || is.logical(data))
    stopifnot(is.numeric(start) || is.logical(start))
    stopifnot(is.numeric(count) || is.logical(count))
    stopifnot(is.logical(pack))

    stopifnot(isTRUE(na.mode %in% c(0,1,2))) 

    #-- Inquire the variable ---------------------------------------------------#
    varinfo <- try(var.inq.nc(ncfile, variable))

    if(class(varinfo) == "try-error" || is.null(varinfo))
        return(invisible(NULL))

    ndims <- varinfo$ndims

    if ((is.character(data) || is.raw(data)) && varinfo$type != "NC_CHAR") {
        stop("R character data can only be written to NC_CHAR variable",call.=FALSE)
    }

    #-- Get the varid as integer if necessary ----------------------------------#
    ifelse(is.character(variable), varid <- varinfo$id, varid <- variable)

    #-- Get correct mode (numeric/character) if data contain only NAs ----------#
    if(is.logical(data)) {
	if(varinfo$type == "NC_CHAR")
	    mode(data) <- "character"
	else
	    mode(data) <- "numeric"
    }

    #-- Check length of character strings --------------------------------------#
    if (is.character(data)) {
      if (ndims > 0) {
        strlen <- dim.inq.nc(ncfile, varinfo$dimids[1])$length
      } else {
        strlen <- 1
      }
      if (max(nchar(data,type="bytes")) > strlen) {
        stop("String length exceeds netcdf dimension", call.=FALSE)
      }
    }

    #-- Replace NA by dimensions of data ---------------------------------------#
    if (any(is.na(start))) {
      start <- rep(1, ndims)
    }

    if (any(is.na(count))) {
      if (!is.null(dim(data))) {
	count <- dim(data)
      } else if (ndims > 0) {
	count <- length(data)
      } else {
	count <- integer(0)
      }
      if (is.character(data) && ndims > 0) {
        count <- c(strlen, count)
      }
    }

    if(length(start) != ndims || length(count) != ndims) {
      stop("Length of start/count is not ndims", call.=FALSE)
    }

    #-- Check that length of data and count match ------------------------------#
    if (is.character(data) && ndims > 0) {
      numelem <- prod(count[-1])
    } else {
      numelem <- prod(count)
    }
    if (length(data) != numelem) {
      stop("Mismatch between count and length(data)", call.=FALSE)
    }

    #-- Pack variables if requested (missing values are preserved) -------------#
    if(pack && is.numeric(data)) {
        offset <- try(att.inq.nc(ncfile, varinfo$name, "add_offset"),
		      silent=TRUE)
        scale  <- try(att.inq.nc(ncfile, varinfo$name, "scale_factor"),
		      silent=TRUE)
        if((!inherits(offset, "try-error")) &&
	   (!inherits(scale,  "try-error"))) {
	    add_offset   <- att.get.nc(ncfile, varinfo$name, "add_offset")
	    scale_factor <- att.get.nc(ncfile, varinfo$name, "scale_factor")
	    data         <- (data-add_offset)*(1/scale_factor)
        }
    }

    #-- Convert missing value to NA if defined in NetCDF file ------------------#
    if(is.numeric(data) && any(is.na(data))) {
      na.flag <- 0

      missval.flag <- 0
      fillval.flag <- 0

      fillval <- try(att.inq.nc(ncfile, varinfo$name, "_FillValue"   ), 
	  silent=TRUE)
      missval <- try(att.inq.nc(ncfile, varinfo$name, "missing_value"),
	  silent=TRUE)

      if(!(class(fillval) == "try-error"))
	  if(!is.null(fillval))
	      fillval.flag <- 1
      if(!(class(missval) == "try-error"))
	  if(!is.null(missval))
	      missval.flag <- 1

      if(na.mode == 0 && missval.flag == 1) {
	  na.value <- att.get.nc(ncfile, varinfo$name, "missing_value")
	  na.flag  <- 1
      }
      if(na.mode == 0 && fillval.flag == 1) {
	  na.value <- att.get.nc(ncfile, varinfo$name, "_FillValue")
	  na.flag  <- 1
      }

      if(na.mode == 1 && fillval.flag == 1) {
	  na.value <- att.get.nc(ncfile, varinfo$name, "_FillValue")
	  na.flag  <- 1
      }

      if(na.mode == 2 && missval.flag == 1) {
	  na.value <- att.get.nc(ncfile, varinfo$name, "missing_value")
	  na.flag  <- 1
      }

      if(na.flag == 1)
	  data[is.na(data)] <- as.numeric(na.value)
      else
	  stop("Found NAs but no missing value attribute", call.=FALSE)
    }

    #-- Switch from R to C convention ------------------------------------------#
    c.start <- rev(start-1)
    c.count <- rev(count)

    #-- C function calls -------------------------------------------------------#
    if(is.numeric(data)) {
        if (!is.double(data)) {
          data <- as.double(data)
        }
	nc <- .Call("R_nc_put_vara_double",
        	    as.integer(ncfile),
		    as.integer(varid),
		    as.integer(c.start),
		    as.integer(c.count),
		    as.integer(ndims),
		    data,
		    PACKAGE="RNetCDF")
    } else {
        stopifnot(is.character(data) || is.raw(data))
        nc <- .Call("R_nc_put_vara_text",
        	    as.integer(ncfile),
		    as.integer(varid),
		    as.integer(c.start),
		    as.integer(c.count),
		    as.integer(ndims),
                    as.integer(is.raw(data)),
		    data,
		    PACKAGE="RNetCDF")
    }
		    
    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  var.rename.nc()                                                              #
#-------------------------------------------------------------------------------#

var.rename.nc <- function(ncfile, variable, newname)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF") 
    stopifnot(is.character(variable) || is.numeric(variable))
    stopifnot(is.character(newname))
    
    #-- Look if handle variable by name or ID ----------------------------------#
    varid   <- -1
    varname <- ""

    ifelse(is.character(variable), nameflag <- 1, nameflag <- 0)
    ifelse(is.character(variable), varname <- variable, varid <- variable)

    #-- C function call --------------------------------------------------------#
    nc <- .Call("R_nc_rename_var",
        	as.integer(ncfile),
		as.integer(varid),
		as.character(varname),
		as.integer(nameflag),
		as.character(newname),
		PACKAGE="RNetCDF")

    if(nc$status != 0)
        stop(nc$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  read.nc()                                                                    #
#-------------------------------------------------------------------------------#

read.nc <- function(ncfile, unpack=TRUE)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(class(ncfile) == "NetCDF")
    stopifnot(is.logical(unpack))

    #-- Initialise storage -----------------------------------------------------#
    nvars    <- file.inq.nc(ncfile)$nvars
    varnames <- character(nvars)
    retlist  <- vector("list", nvars)

    #-- Read data from each variable -------------------------------------------#
    for(i in seq_len(nvars)) {
        retlist[[i]] <- var.get.nc(ncfile, i-1, unpack=unpack)
        varnames[i]  <- var.inq.nc(ncfile, i-1)$name
    }

    #-- Set names of list elements ---------------------------------------------#
    names(retlist) <- varnames

    return(retlist)
}


#===============================================================================#
#  Udunits library functions							#
#===============================================================================#

#-------------------------------------------------------------------------------#
#  utcal.nc()                                                                   #
#-------------------------------------------------------------------------------#

utcal.nc <- function(unitstring, value, type="n")
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(is.character(unitstring))
    stopifnot(is.numeric(value) && !any(is.na(value)))
    stopifnot(type == "n" || type =="s" || type == "c" )
    
    count <- length(value)
   
    #-- C function call to udunits calendar function -----------------------#
    ut <- .Call("R_ut_calendar", 
	        as.character(unitstring), 
		as.integer(count),
		as.double(value),
		PACKAGE="RNetCDF")

    #-- Return object if no error ------------------------------------------#
    if(ut$status == 0) {
	if(type == "n") {
	    colnames(ut$value) <- c("year", "month", "day", "hour", 
		"minute", "second")
	    return(ut$value)
	} else if (type == "s") {
	    x <- apply(ut$value, 1, function(x){paste(x[1],"-",
		       sprintf("%02g",x[2]),"-",sprintf("%02g",x[3])," ",
		       sprintf("%02g",x[4]),":",sprintf("%02g",x[5]),":",
		       sprintf("%02g",x[6]),sep="")})
	    return(x)
        } else if (type == "c") {
            ct <- as.POSIXct(
                        utinvcal.nc("seconds since 1970-01-01 00:00:00 +00:00",ut$value),
                        tz="UTC", origin=ISOdatetime(1970,1,1,0,0,0,tz="UTC"))
            return(ct)
        }
    } else {
	stop(ut$errmsg, call.=FALSE)
    }
}


#-------------------------------------------------------------------------------#
#  utinit.nc()                                                                  #
#-------------------------------------------------------------------------------#

utinit.nc <- function(path="")
{
    ut <- .Call("R_ut_init", 
                as.character(path),
		PACKAGE="RNetCDF")
		
    if(ut$status != 0)
        stop(ut$errmsg, call.=FALSE)
}


#-------------------------------------------------------------------------------#
#  utinvcal.nc()                                                                #
#-------------------------------------------------------------------------------#

utinvcal.nc <- function(unitstring, value)
{
    #-- Check args -------------------------------------------------------------#
    stopifnot(is.character(unitstring))

    if (is.character(value)) {
	stopifnot(all(nchar(value) == 19))
	value <- cbind(substr(value,1,4),
		       substr(value,6,7),
		       substr(value,9,10),
		       substr(value,12,13),
		       substr(value,15,16),
		       substr(value,18,19))

	value <- matrix(as.numeric(value),ncol=6)
    } else if (inherits(value,"POSIXct")) {
        value <- utcal.nc("seconds since 1970-01-01 00:00:00 +00:00",
                     as.numeric(value), 'n')
    }

    stopifnot(is.numeric(value) && !any(is.na(value)))

    count <- length(value)
    
    if(is.vector(value) && count %% 6  != 0)
	stop("length(value) not divisible by 6", call.=FALSE)

    if(is.matrix(value) && ncol(value) != 6) 
	stop("ncol(value) not 6", call.=FALSE)
    
    #-- C function call --------------------------------------------------------#
    ut <- .Call("R_ut_inv_calendar", 
		as.character(unitstring), 
		as.integer(count),
		as.double(value),
		PACKAGE="RNetCDF")

    #-- Return object if no error ----------------------------------------------#
    if(ut$status == 0) {
	return(ut$value)
    } else {
	stop(ut$errmsg, call.=FALSE)
    }
}


#===============================================================================#

#===============================================================================#
#  SCRATCH									#
#===============================================================================#

