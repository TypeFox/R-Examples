
#===============================================================
ncdim_same <- function( d1, d2 ) {

	if( class(d1) != "ncdim4" ) 
		stop("error, class of first passed argument is not ncdim!")
	if( class(d2) != "ncdim4" ) 
		stop("error, class of first passed argument is not ncdim!")

	if( d1$name != d2$name )
		return(FALSE)

	if( d1$len != d2$len )
		return(FALSE)

	if( d1$unlim != d2$unlim )
		return(FALSE)

	return(TRUE)
}

#====================================================================================================
# This is the private interface that actually does the 
# netCDF calls.  User code should never go through this.
# To make a ncdim object, use ncdim_def() instead.
# This makes BOTH the dim AND the dimvar (and RETURNS
# dimid AND dimvarid).
#
ncdim_create <- function( nc, d, verbose=FALSE ) {  

	if( class(nc) != "ncdf4" ) 
		stop("ncdim_create: passed nc NOT of class ncdf4!")
	if( verbose )
		print(paste("ncdim_create: entering for ncid=",nc$id))

	if( class(d) != "ncdim4" ) 
		stop("ncdim_create: passed d NOT of class ncdim4!")
	if( verbose )
		print(paste("ncdim_create: entering for dim",d$name))

	#-----------------------------------------------------------------------
	# Figure out the ncid to use.  If this dim is in the root group, it will
	# just be the ncid.  Otherwise, it will be the group id
	#-----------------------------------------------------------------------
	if( nslashes_ncdf4( d$name ) == 0 )
		gidx <- 1
	else
		{
		dims_fqgn <- nc4_basename( d$name, dir=TRUE )
		gidx      <- nc$fqgn2Rindex[[ dims_fqgn ]]
		if( is.null(gidx))
			stop(paste("internal error: did not find dim's fully qualified group name", dims_fqgn," in list of groups for file", nc$filename))
		}
	ncid2use  <- nc$group[[gidx]]$id

	#-------------------
	# Make the dimension
	#-------------------
	ncdim <- list()
	ncdim$error <- -1
	ncdim$id    <- -1	# this is only briefly an integer, changed to ncid object below
	sizetouse   <- d$len
	name2use    <- nc4_basename( d$name ) 
	if( d$unlim )
		sizetouse <- 0
	if( verbose )
		print(paste("ncdim_create: about to call R_nc4_def_dim with ncid2use=", 
			ncid2use, " NOT fully qualified dim name=",name2use))
	ncdim<-.C("R_nc4_def_dim",
		as.integer(ncid2use),
		as.character(name2use),
		as.integer(sizetouse),
		id=as.integer(ncdim$id),
		error=as.integer(ncdim$error),
		PACKAGE="pbdNCDF4")
	if( ncdim$error != 0 ) 
		stop("Error in dim.create.ncdf!")
	#------------------------------------------------------
	# NOTE that ncdim$id is just an ORDINARY INTEGER, not a
	# ncid object!  The 'ncdim' structure is not returned,
	# it's just used locally in this routine.
	#------------------------------------------------------

	#-----------------------------
	# Make the dimvar if requested
	#-----------------------------
	dimvar<-list()
	if( d$create_dimvar ) {
		if( verbose ) print(paste("ncdim_create: making dimvar for dim",d$name))
		dimvar$id    <- -1	# this is only briefly an integer, changed to ncid object below
		dimvar$error <- -1
		if( storage.mode(d$vals) == "integer" ) {
			if( verbose )
				print(paste("ncdim_create: about to call R_nc4_def_var_int for dimvar",name2use))
			dimvar<-.C("R_nc4_def_var_int",
				as.integer(ncid2use),
				as.character(name2use),
				as.integer(c(1)),
				as.integer(ncdim$id),	
				id=as.integer(dimvar$id),
				error=as.integer(dimvar$error),
				PACKAGE="pbdNCDF4")
			}
		else
			{
			if( verbose )
				print(paste("ncdim_create: about to call R_nc4_def_var_double for dimvar",name2use))
			dimvar<-.C("R_nc4_def_var_double",
				as.integer(ncid2use),
				as.character(name2use),
				as.integer(c(1)),
				as.integer(ncdim$id), 
				id=as.integer(dimvar$id),
				error=as.integer(dimvar$error),
				PACKAGE="pbdNCDF4")
			}
		if( dimvar$error != 0 ) 
			stop("Error defining dimvar in routine dim.create.ncdf")
		#-------------------------------------------------------
		# NOTE that dimvar$id is just an ORDINARY INTEGER, not a
		# ncid object!  The 'dimvar' structure is not returned,
		# it's just used locally in this routine.
		#-------------------------------------------------------

		#---------------------------------
		# Put in the dimvals as specified.
		#---------------------------------
		nc_enddef( nc )		# Must exit define mode for this
		rv <- list()
		rv$error <- -1
		start <- 0		# Use C convention
		count <- length(d$vals)
		if( count > 0 ) {
			if( storage.mode(d$vals) == "integer" ) {
				if( verbose )
					print(paste("ncdim_create: about to call R_nc4_put_vara_int dimvals for dimvar",d$name))
				rv <- .C("R_nc4_put_vara_int",
					as.integer(ncid2use),
					as.integer(dimvar$id),
					as.integer(start),
					as.integer(count),
					as.integer(d$vals),
					error=as.integer(rv$error),
					PACKAGE="pbdNCDF4")
				}
			else if( storage.mode(d$vals) == "double" ) {
				if( verbose )
					print(paste("ncdim_create: about to call R_nc4_put_vara_double dimvals for dimvar",d$name))
				rv <- .C("R_nc4_put_vara_double",
					as.integer(ncid2use),
					as.integer(dimvar$id),
					as.integer(start),
					as.integer(count),
					as.double(d$vals),
					error=as.integer(rv$error),
					PACKAGE="pbdNCDF4")
				}
			else
				stop(paste("ncdim_create: unknown storage mode:",storage.mode(d$vals),"for dim",d$name))
			if( rv$error != 0 )
				stop("Error in ncdim_create, while writing dimvar values!")
			}
		nc_redef( nc )	# Go back into define mode

		#----------------------------------------------------
		# Set the dimension's (dimvar's, actually) attributes
		#----------------------------------------------------
		if( (! is.null(d$units)) && (nchar(d$units)>0))
			ncatt_put_inner( ncid2use, dimvar$id, "units",     d$units,    definemode=TRUE, verbose=verbose )
		if( (! is.null(d$longname)) && (nchar(d$longname)>0))
			ncatt_put_inner( ncid2use, dimvar$id, "long_name", d$longname, definemode=TRUE, verbose=verbose )

		if( (! is.null( d$calendar )) && (nchar(d$calendar)>0))
			ncatt_put_inner( ncid2use, dimvar$id, "calendar", d$calendar, definemode=TRUE, verbose=verbose )

		}	# end of "if(create_dimvar)"
	else
		{
		if( verbose ) print(paste("ncdim_create: NOT making dimvar for dim",d$name))
		#----------------------------------------------------------
		# if we were NOT asked to create the dimvar (via an empty
		# units string) than make sure NO dim values were specified
		# except simple integers from 1 to len!
		#----------------------------------------------------------
		if( (storage.mode( d$vals ) != "integer" ) || (d$vals[1] != 1) || (d$vals[d$len] != d$len))
			stop(paste("Error trying to create dimension named",d$name,": the passed units string",
				"was empty, which indicates that NO dimensional variable is to be created.",
				"In this case, the dimension values MUST be simple integers from 1 to the length",
				"of the dimension (e.g., 1:len)"))
		dimvar$id = -1
		}

	if( verbose )
		print(paste("ncdim_create: exiting for ncid=",nc$id,"  dim=",d$name))

	#----------------------------------------------
	# Return the dimvar ID of the newly created dim
	#----------------------------------------------
	return(c(ncdim$id,dimvar$id,ncid2use,gidx))
}

#===============================================================
# NOTE that this does NOT return a full-fledged "ncdim" object,
# it is the lower-level interface and only returns the portions
# of a ncdim object that directly correspond to a real netCDF
# dimension.  The ncdim object, by contrast, also has information
# from the dimvar.
#
# Internal use only.
#
ncdim_inq <- function( ncid, dimid ) {

	#------------------------------------------------------
	# This is a string NC_MAX_LEN long, to provide storage.
	# I'm not sure if this is needed or not....
	#------------------------------------------------------
	str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"

	rv <- list()
	rv$dimname <- str.nc.max.name
	rv$dimlen  <- -1
	rv$error   <- -1
	rv$unlim   <- 0
	rv <- .C("R_nc4_inq_dim",
		as.integer(ncid),
		as.integer(dimid),
		dimname=as.character(rv$dimname),
		dimlen=as.integer(rv$dimlen),
		unlim=as.integer(rv$unlim),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) {
		stop(paste("using ncid ",ncid," dimid ",dimid))
		}
	d <- list()
	d$name  <- rv$dimname
	d$len   <- rv$dimlen
	d$unlim <- (rv$unlim == 1)

	return(d)
}

#===============================================================
# Internal use only
#
# Returns -1 if the dim is NOT found in the file or group, and the
# dimid of the dim otherwise.
#
ncdim_id <- function( nc, dimname ) {

	if( mode(nc) != 'numeric' )
		stop("error, must be passed a numeric first arg: ncid2use")

	if( mode(dimname) != 'character' )
		stop("Error, must be passed a character second arg: dimname" )

	rv       <- list()
	rv$dimid <- -1
	rv <- .C("R_nc4_inq_dimid", 
		as.integer(nc),
		as.character(dimname),
		dimid=as.integer(rv$dimid),
		PACKAGE="pbdNCDF4")
	if( rv$dimid != -1 )
		rv$dimid <- rv$dimid 
	return(rv$dimid)
}

