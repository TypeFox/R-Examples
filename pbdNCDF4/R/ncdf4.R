# This provides an interface to netCDF-4 functions.  As a general
# comment, awkwardness arises from two things:
#	1. R starts counting at 1, and netCDF counting
#	   starts at 0.
#	2. R array subscripts go in Fortran order (XYZT),
#	   while netCDF subscripts go in C order (TZYX).
# We take care of these problems EXCLUSIVELY in the R code,
# NOT in the C interface code!!!   This means that it is
# always the responsibility of the R program to take these
# differences into account.  From the point of view of an
# R program that calls any of these functions, they are 
# strictly R compliant (Fortran order, counting starts at 1).
#
# David W. Pierce
# Climate Research Division
# Scripps Institution of Oceanography
# dpierce@ucsd.edu
# 15-January-2010
#
#-----------------------------------------------------------------
# Netcdf-version 4 notes.
# Not much needs to be changed in the interface to make it work
# for netcdf version 4.  Basically, if you have a slash in the
# variable or dim name, such as "models/Temperature", then "models"
# is a group that will be auto-created as necessary.  Same goes
# for dims.  NOTE that this is fundamentally a different way
# of doing things than the raw netcdf C interface does.  The
# C interface uses DIFFERENT ncid's, corresponding to each group.
# THIS interface always uses the ncdf object, but the name of
# the dim or var must be fully qualified if there is an 
# ambiguity about which one it refers to.
#
# Once change I did make was that becuase ID's are now more of
# an abstract concept, and ID numbering can be arbitrary, the
# USER interface will ONLY take the netcdf objects, not integers.
# For example, varid's used to be numbered, in C convention, 
# from 0 to N-1.  Now they can be arbitrary.  As a result,
# all the "ids" stored in the ncdf object are the actual C
# values themselves.  We no longer need to, or try to, adjust
# between "R" values and "C" values.
#
# varid's have been changed to reflect this.  Now, varids are
# an object of type 'ncid'.  That object holds the "raw" varid, 
# plus the netcdf group id needed to access the variable.  
#
# dimid's, on the other hand, are still simple integers.  This
# seems to be a characteristic of netcdf files, possibly because
# dims can be seen if they exist in the group OR ANY PARENT GROUP.
# Obviously enforcing that would be a headache, which is made
# easier if you simply have dims numbered sequentially across
# the whole file (as opposed to vars, which are numbered 
# sequentially, starting at 0, in EACH GROUP).
#
# In netcdf-4, there can be multiple unlimited dims.  In order to
# avoid breaking existing code, nc$unlimdimid is set to the FIRST
# unlimited dim in the file.  Technically, unlimdimid is not set
# to an ID, it's set to an index into the global nc$dim[[]] list.
#-----------------------------------------------------------------
#
# Here are the relevant objects:
#
# class: ncdf4 is a list with the following fields:
#	filename: name of the file, or "IN-MEMORY"
#	id	: netcdf file id
#	ndims	: integer # of dims in the file
#	nvars	: integer # of vars in the file that are **NOT** dimvars
#	natts	: integer # of global attributes
#	group   : list of groups in the file.  group[[1]] is always the root group
#	ngroups : integer # of groups; will always be at least 1, indicating the root group;
#		this is just length(group)
#	dim	: a list of ncdim objects
#	var	: a list of ncvar objects
#	writable: TRUE or FALSE
#
# class: ncdim4 (returned by ncdim_def, which creates a NEW 
#		netCDF dimension in memory, and part of the list of dims
#		in a ncdf object. NOTE that this is NOT what is returned
#		by ncdim_inq, which is the low-level netCDF dim, not
#		the user-level R version of a netCDF dim.
#     *	name	: character dim name
#	units	: character units in udunits format.  
#	vals	: a vector of dimension values
#     *	len	: size of this dimension
#	calendar : if the dimvar has a 'calendar' attribute, this holds its value
#	id	: NEW IN VERION 4: this is an opaque object of type 'ncdimid'
#	dimvarid: IFF this dim corresponds to a netCDF dim in an existing file,
#		  then this field will hold the dimvarid, or -1 if no dimvar.
#     *	unlim	: boolean, T or F to indicated unlimited or not
#	create_dimvar : usually TRUE; if FALSE, then no dimvar will be created,
#		  AND the units string must be empty, AND the values must
#		  be simple integers from 1 to the length of the dim.
#	longname : the 'long_name' attribute for the DIMVAR, or defaults to
#		  name if there is no dimvar, OR the long_name attribute 
#		  is not set.
#
# *=Indicates element is filled out by the low-level routine 'ncdim_inq'.
#   Other elements are filled out in 'open.ncdf'.
#
#
# class: ncvar4 (returned by ncvar_def, which creates a NEW
#		netCDF variable in memory, and part of the list of vars
#		in a ncdf object.  NOTE that this is NOT what is
#		returned by ncvar_inq, which is the low-level netCDF
#		var, NOT the user-level R version of a netCDF var.
#	name	: character var name
#	units	: character units in udunits format
#	missval	: the 'missing_value' attribute, or defaults to default_missval_ncdf4(). NOTE: 
#		  if the var has no missing value then this is 'NA'.  This does not mean
#		  that the missing value is NA -- the missing value cannot be NA.  It means
#		  the variable has no missing value.  For instance, character variables have
#		  no default missing value.
#	longname: the 'long_name' attribute, or defaults to name
#	id	: the varid of this variable, IFF it is in a file already
#	ndims	: number of dims this variable has
#	dim	: a list of type ncdim, which is this variable's dims
#	unlim	: boolean, T if this var has an unlimited dim, F otherwise
#	varsize : a convenience array that gives the (X,Y,Z,T) size of
#		  the variable.
#	prec    : The precision of the ON-DISK representation of the
#		  variable.  Can be "byte", "short", "float", "double", 
#	    	  "integer", or "char".
#	natts	: integer # of this variable's attributes
#
#======================================================================================================
nc_version <- function() {
	
	return("ncdf4_1.8_20120922")

}

#====================================================================================================
# This is the public interface for making a netCDF dimension
# class.  It makes it in memory, not in a file!  Returns a
# "ncdim" object.  The dimvar will have the same precision
# as the passed values.  Therefore, it will always have double
# precision unless the vals are passed like this: as.integer(vals),
# in which case they will be integer.
#
# Example useage:
#	lon <-  ncdim_def("Lon", "degreesE", 0:359)
#	lat <-  ncdim_def("Lat", "degreesN", -90:90)
#	time <- ncdim_def("time", "days since 1900-01-01", 0, unlim=T)
#
# Netcdf-4 useage:
#	lon <- ncdim_def( "models/pcm/Lon", ... )
# where groups are indicated by things before the slashes.
# 
# NOTE that fully qualified dim names NEVER start with a slash!  This
# is required for backward compatability with old code.
#
ncdim_def <- function( name, units, vals, unlim=FALSE, create_dimvar=TRUE, 
		calendar=NA, longname=name ) {

	#---------------
	# Validate input
	#---------------
	if( (! is.character(name)) || (nchar(name) < 1))
		stop("Passed a dim name that is NOT a string of characters!")

	if( ! is.character(units) ) 
		stop("Passed a dim units that is NOT a string of characters!")

	if( ! is.numeric(vals))
		stop("Only numeric dimension values are supported")

	if( ! is.na(calendar)) {
		if( (! is.character(calendar)) || (nchar(calendar)<1))
			stop("Passed a calendar that is NOT a string of characters!")
		}

	if( substr(name,1,1) == '/' )
		stop(paste("Error: passed dim name", name, "starts with a slash.  Fully qualified dim names can NEVER start with a slash (this is required for backwards compatability).  Leave off any leading slash!"))

	len <- length(vals)
	if( ! create_dimvar ) {
		if( (units != '') || (storage.mode(vals) != "integer" ) || (vals[1] != 1) || (vals[len] != len))
			stop(paste("Error trying to create dimension named",name,": create_dimvar was specified",
				"to be FALSE, which indicates that NO dimensional variable is to be created;",
				"in this case, the unit string MUST be empty ('') and the dimension values MUST",
				"be simple integers from 1 to the length of the dimension (e.g., 1:len)"))
		}

	dim <- list()
	dim$name 	  <- name		# Note this is a fully qualified dim name and NEVER starts with a slash
	dim$longname   	  <- longname
	dim$units  	  <- units
	dim$vals   	  <- vals
	dim$len    	  <- len

	#---------------------------------------------------------------------------
	# In netcdf v4, there can be multiple vars with the same name and same ID,
	# if they are in different groups.  So we have to have a more complex object
	# for the dimvar ID than just a simple integer.  This does not seem to be
	# true for dims, which still appear to have unique ID's across all groups
	# in the file.
	#
	# Elements of the dimvar's 'ncid' object:
	#	id: the integer to use to access this dim when calling the 
	#		C routines to access the dim WITH THE CORRECT GROUP ID
	#	group_id: the group ID to use when accessing the C routines to
	#		manipulate this dim
	#	group_index: the index in to the nc$group list for the group where
	#		this dim lives
	#	list_index: index into the global nc$var[[]] list for this var.
	#		NOTE that dimvars are not on this list.
	#	isdimvar: is TRUE if the var is a dimvar, FALSE otherwise.
	# Note that none of these is valid until the dim is actually created on disk
	#---------------------------------------------------------------------------
	dim$id 		  <- -1
	dim$dimvarid 	  <- ncdf4_make_id( isdimvar=TRUE )	# make an EMPTY ncid object

	dim$unlim  	  <- unlim
	if( !is.na(calendar)) {
		if( ! create_dimvar ) {
			print(paste("Error in ncdim_def while trying to define dimension",name,":"))
			print(paste("This dimension was specified to have a calendar (",calendar,") BUT"))
			print(paste("also it was specified that no dimvar would be created, via create_dimvar=FALSE!"))
			print(paste("This is a contradiction since the calendar attribute must be stored"))
			print(paste("on the dimvar.  Either get rid of the calendar, or allow a dimvar"))
			print(paste("to be created."))
			stop("Halting.")
			}
		dim$calendar <- calendar
		}
	dim$create_dimvar <- create_dimvar
	attr(dim,"class") <- "ncdim4"

	return(dim)
}

#====================================================================================================
# Produce a more useful listing of the netcdf file than
# just a dump of the object.
#
print.ncdf4 <- function( x, ... ) {

	nc <- x

	is_netcdf_v4 = (nc$format == 'NC_FORMAT_NETCDF4')

	print(paste("File ",nc$filename, " (", nc$format, "):", sep=''))
	print("")
	print(paste("    ",nc$nvars,"variables (excluding dimension variables):"))
	if( nc$nvars > 0 ) {
		for( i in 1:nc$nvars ) {
			nd <- nc$var[[i]]$ndims
			dimstring <- '['
			if( nd > 0 ) {
				for( j in 1:nd ) {
					dimstring <- paste(dimstring,nc$var[[i]]$dim[[j]]$name,sep='')
					if( j < nd )
						dimstring <- paste(dimstring,',',sep='')
					}
				}
			dimstring <- paste(dimstring,'] ',sep='')

			chunk_tag = ''
			compress_tag = ''
			if( is_netcdf_v4 ) {

				#----------------------------
				# Handle chunking information
				#----------------------------
				if( is.null(nc$var[[i]]$storage) || nc$var[[i]]$storage == 1 )
					chunk_tag = "  (Contiguous storage)"
				else
					{
					chunk_tag = "  (Chunking: ["
					for( j in 1:nd ) {
						chunk_tag = paste( chunk_tag, nc$var[[i]]$chunksizes[j], sep='' )
						if( j < nd )
							chunk_tag = paste( chunk_tag, ",", sep='' )
						}
					chunk_tag = paste( chunk_tag, "])", sep='' )
					}

				#---------------------------------------
				# Handle shuffle/compression information
				#---------------------------------------
				is_shuffle  = (nc$var[[i]]$shuffle == 1)
				is_compress = (!is.na(nc$var[[i]]$compression))
				if( (!is_shuffle) && (!is_compress))  
					compress_tag = ""
				else if( is_shuffle && (!is_compress))
					compress_tag = "(Compression: shuffle)"
				else if( (!is_shuffle) && is_compress )
					compress_tag = paste("(Compression: level ", nc$var[[i]]$compression, ")", sep='' )
				else
					compress_tag = paste("(Compression: shuffle,level ", nc$var[[i]]$compression, ")", sep='' )
				}
			print(paste("        ", nc$var[[i]]$prec, ' ', nc$var[[i]]$name, dimstring, chunk_tag, "  ", compress_tag, sep='' ))
			atts <- ncatt_get( nc, nc$var[[i]]$name )
			natts <- length(atts)
			if( natts > 0 ) {
				nms <- names( atts )
				for( ia in 1:natts ) 
					print(paste("            ", nms[ia], ": ", atts[[ia]], sep='' ))
				}
			}
		}

	print("")
	print(paste("    ",nc$ndims,"dimensions:"))
	if( nc$ndims > 0 )
		for( i in 1:nc$ndims ) {
			tag <- ''
			if( nc$dim[[i]]$unlim )
				tag <- '   *** is unlimited ***'
			print(paste("        ", nc$dim[[i]]$name, "  Size:", nc$dim[[i]]$len, tag, sep='' ))
			atts <- ncatt_get( nc, nc$dim[[i]]$name )
			natts <- length(atts)
			if( natts > 0 ) {
				nms <- names( atts )
				for( ia in 1:natts ) 
					print(paste("            ", nms[ia], ": ", atts[[ia]], sep='' ))
				}
			}

	#--------------------------
	# Now get global attributes
	#--------------------------
	atts <- ncatt_get( nc, 0 )
	natts <- length(atts)
	if( natts > 0 ) {
		print('')
		print(paste('    ', natts, ' global attributes:', sep=''))
		nms <- names( atts )
		for( ia in 1:natts ) 
			print(paste("        ", nms[ia], ": ", atts[[ia]], sep='' ))
		}
}

#==========================================================================================================
# This is the public interface for making a netCDF variable
# class.  It makes it in memory, not in a file!  Returns
# a "ncvar" object.
#
# Example useage, where "lon" and "lat" are objects of class "ncdim":
#
#	ncvar <- ncvar_def( "temp", "degC", list(lon,lat), 1.e30, "Temperature" )
# or
#	ncvar <- ncvar_def( "temp", "degC", lon, 1.e30, "Temperature" )
#
# Example useage for netcdf library version 4:
#
#	ncvar <- ncvar_def( "models/pcm/Temperature", ... )
#
# NOTE that the passed dimentions (in this example, "lon" and "lat"
# should have been made by "ncdim_def", and so are of class
# "ncdim".   Argument "prec" will indicate what precision the variable
# is made on the disk.  Allowed values are "short", "integer", "float", 
# "double", "byte", and "char".
# 'missval' is the value to be assigned to the 'missing_value' attribute
# for the variable.  It should be a number representable in the precision-type
# of the variable.  So, for example, for single precision, a floating point
# number with magnitude less than 1.e36 should be used.  As a special value,
# if missval is set to NA, then *no* missing value will be created for the 
# variable.
#
# To make a var with no dims, pass an empty list: "list()"

ncvar_def <- function( name, units, dim, missval, longname=name, prec="float", 
		shuffle=FALSE, compression=NA, chunksizes=NA, verbose=FALSE ) {

	if( verbose ) print('ncvar_def: entering')

	#-----------------------------------
	# Check inputs for reasonable values
	#-----------------------------------
	if( (! is.character(name)) || (nchar(name) < 1))
		stop("Passed a var name that is NOT a string of characters!")

	if( ! is.character(units))	# Note: an empty string is OK 
		stop("Passed a var units that is NOT a string of characters!")

	if( (! is.character(prec)) || (nchar(prec) < 1))
		stop("Passed a var precision (prec) that is NOT a string of characters!")
	if( prec == 'single' )
		prec = 'float'

	if( verbose ) print(paste('ncvar_def: prec=', prec))
	if( !is.na(compression)) {
		if(!is.numeric(compression))
			stop("Compression parameter, if supplied, must be an integer between 1 (least compression) and 9 (most compression)") 
		compression = as.integer(compression)
		if((compression < 1) || (compression > 9))
			stop("Compression parameter, if supplied, must be an integer between 1 (least compression) and 9 (most compression)") 
		}

	if( (length(chunksizes)>1) || (!is.na(chunksizes))) {
		if( (! is.numeric(chunksizes)) || (length(chunksizes) != length(dim)))
			stop("Chunksizes parameter, if supplied, must be a vector of integers with length equal to the dimension of the variable")
		}

	if( shuffle ) {
		if( (prec == "float") || (prec == "double") || (prec == "single") )
			print(paste("Warning: shuffle is turned on for variable", name, "but that var is of precision", prec,
				"and shuffle ONLY has an effect for integer variables."))
		if( is.na(compression))
			print(paste("Warning: shuffle is turned on for variable", name, "but compression is NOT turned on. Quoting from the netCDF docs: 'There is no benefit from using the shuffle filter without also using compression.'"))
		}	

	if( substr(name,1,1) == '/' )
		stop(paste("Error: passed var name", name, "starts with a slash.  Fully qualified var names can NEVER start with a slash (this is required for backwards compatability).  Leave off any leading slash!"))

	#----------------------
	# Make our ncvar object
	#----------------------
	if( verbose ) print(paste('ncvar_def: making ncvar object for var', name ))
	var <- list()
	var$name     	<- name 	# Note this is a fully qualifed var name and never starts with a slash 
	var$units    	<- units
	var$longname 	<- longname
	var$shuffle  	<- shuffle
	var$compression <- compression
	var$chunksizes  <- chunksizes	# a vector of integers; the length of the vector == ndims

	#----------------------------------------------------------------------------------
	# New: in 1.1, you can specify NOT to have a missing value by NOT passing a missing
	# value argument at all. Also, you can now specify a missing value of NA by passing
	# a NA.
	#----------------------------------------------------------------------------------
	if( missing(missval) || is.null(missval))
		var$make_missing_value = FALSE
	else
		{
		var$make_missing_value = TRUE
		var$missval  	<- missval
		if( storage.mode(missval) == "character" )
			prec <- 'char'
		}

	#---------------------------------------------------------------------------
	# In netcdf v4, there can be multiple vars with the same name and same ID,
	# if they are in different groups.  So we have to have a more complex object
	# for the "ID" than just a simple integer.
	# Elements:
	#	id: the integer to use routines to access this vars when calling the 
	#		C routines to access the dim WITH THE CORRECT GROUP ID
	#	group_id: the group ID to use when accessing the C routines to
	#		manipulate this vars
	#	group_index: the index in to the nc$group list for the group where
	#		this vars lives
	#	list_index: index into the glocal nc$var[[]] list that corresponds
	#		to this var
	#	isdimvar: TRUE if this is a dimvar, FALSE otherwise
	# Note that none of these is valid until the vars is actually created on disk
	#---------------------------------------------------------------------------
	var$id 		  <- ncdf4_make_id()	# make an EMPTY ncid object

	attr(var,"class") <- "ncvar4"

	if( (prec != "short")   && (prec != "float") && (prec != "double") && 
	    (prec != "integer") && (prec != "char")   && (prec != "byte"))
		stop(paste("ncvar_def: error: unknown precision specified:",prec,". Known values: short float double integer char byte"))
	var$prec <- prec

	#-----------------------------------------------------
	# Have to figure out if 'dim' is a ncdim object or 
	# a LIST of ncdim objects.
	#-----------------------------------------------------
	if( is.character(class(dim)) && (class(dim) == "ncdim4") )
		dim <- list(dim)
	var$dim  <- dim
	var$ndims <- length(var$dim)
	if( var$ndims > 0 ) {
		for( i in 1:var$ndims ) {
			if( class(var$dim[[i]]) != "ncdim4" ) {
				print(paste("Error, passed variable has a dim that is NOT of class ncdim4!"))
				print(paste("Error occurred when processing dim number",i,"of variable",var$name))
				stop(paste("This dim has class:", class(var$dim[[i]])))
				}
			}
		}

	#-------------------------------------------------
	# A variable is unlimited if any of its dimensions
	# is unlimited
	#-------------------------------------------------
	varunlimited <- FALSE
	if( var$ndims != 0 ) {
		for( i in 1:var$ndims ) {
			if( var$dim[[i]]$unlim ) 
				varunlimited <- TRUE
			}
		}
	var$unlim <- varunlimited

	return(var)
}

#===============================================================================================
# This is the public interface for opening an already-existing
# netCDF file.  If you want to create a NEW netCDF file, use
# nc_create instead.  This returns an object of class 'ncdf',
# used to access netCDF functions.  If you want to modify an 
# already-existing netCDF file, use this function and set write=TRUE.
# By default, this will read in the values of the unlimited dimension(s)
# along with the values of all the other dimensions.  Since this
# can be time comsuming, there is the option of setting 'readunlim' 
# to 'FALSE', which will avoid this behavor.  Then you must read
# in the values yourself with a 'ncvar_get' call, instead of
# accessing the ncid$dim[[DIMNAME]]$vals array.
#
nc_open <- function( filename, write=FALSE, readunlim=TRUE, verbose=FALSE ) {

	if( verbose ) print(paste("nc_open: entering, package version", nc_version() ))

	if( (! is.character(filename)) || (nchar(filename) < 1))
		stop("Passed a filename that is NOT a string of characters!")

	if( verbose )
		print(paste("nc_open: entering, ncdf version=",nc_version()))

	rv <- list()

	if( write )
		rv$cmode <- 1
	else
		rv$cmode <- 0

	rv$id    <- -1
	rv$error <- -1
	rv <- .C("R_nc4_open",
		as.character(filename),
		as.integer(rv$cmode),
		id=as.integer(rv$id),		# note: nc$id is the simple integer ncid of the base file (root group in the file)
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop(paste("Error in nc_open trying to open file",filename))
	if( verbose )
		print(paste("nc_open: back from call to R_nc4_open, ncid=",rv$id))

	#-------------------------------------------------
	# Now we make our elaborate ncdf class object
	#-------------------------------------------------
	nc <- list( filename=filename, writable=write, id=rv$id )
	attr(nc,"class") <- "ncdf4"

	#-------------------------------------------------------------
	# See what format this file is.  Possible (string) values are:
	# 'NC_FORMAT_CLASSIC', 'NC_FORMAT_64BIT', 'NC_FORMAT_NETCDF4',
	# and 'NC_FORMAT_NETCDF4_CLASSIC'
	#-------------------------------------------------------------
	nc$format = ncdf4_format( nc$id )
	if( verbose )
		print(paste("file", filename, "is format", nc$format ))

	#-----------------------------------------------------
	# Get all the groups in the file.  Later we have to 
	# remember that dims and vars can live in groups other
	# than the root group.
	#-----------------------------------------------------
	groups <- list()
	groups[[1]] <- nc_get_grp_info( nc$id, "", nc$format )	# sets $name, $fqgn, $id, $nvars, $ndims, $natts, $ngrps, etc
	if( nc$format == 'NC_FORMAT_NETCDF4' ) {
		gg <- nc_groups_below( groups[[1]], nc$format )
		for( i in nc4_loop(1,length(gg)))
			groups[[1+i]] <- gg[[i]]
		}
	nc$groups <- groups

	if( verbose ) {
		print("Group info:")
		for( ig in 1:length(groups)) {
			print(paste("Group", ig, ": ",
				"name=", groups[[ig]]$name,
				"id=", groups[[ig]]$id,
				"fqgn= \"", groups[[ig]]$fqgn, "\"",
				"nvars=", groups[[ig]]$nvars,
				"ndims=", groups[[ig]]$ndims,
				"dimid=")) 
			print( groups[[ig]]$dimid )
			}
		}

	#---------------------------------------------------------------------------
	# Get general information about the file.  ndims and natts is the sum across
	# all groups in the file.  nvars is set below to only include non-dimvars, 
	# but otherwise is also the sum across all groups.
	#---------------------------------------------------------------------------
	nc$ndims <- 0
	nc$natts <- 0
	tot_nvars_inc_dimvars <- 0	# total number of vars INCLUDING dimvars; regular nc$nvars excludes dimvars
	for( ig in 1:length(groups)) {
		nc$ndims <- nc$ndims + nc$groups[[ig]]$ndims
		nc$natts <- nc$natts + nc$groups[[ig]]$natts
		tot_nvars_inc_dimvars <- tot_nvars_inc_dimvars + nc$groups[[ig]]$nvars
		}

	#--------------------------------------------
	# Get all the dimensions that this file has.
	# Get their values as well (for caching).  
	#--------------------------------------------
	nc$dim        <- list()
	nc$unlimdimid <- -1		# Will be set to the FIRST encountered unlim dim ID
	dimnames      <- character()
	global_dim_counter <- 0		# counter of dims across ALL groups
	for( ig in 1:length(groups)) {

		for( idim in nc4_loop(1,groups[[ig]]$ndims)) {	

			#-----------------------------------------------
			# New in v4: dimids don't go from 0..ndims, they
			# can be arbitrary
			#-----------------------------------------------
			dimid2use <- groups[[ig]]$dimid[idim]
			if( is.na(dimid2use)) {
				print(paste("Error, got a NA as a dimid2use ... group=",
					groups[[ig]]$name," which has ndims=",
					groups[[ig]]$ndims ))
				print("Here are the dimids from the ncgroup object:")
				print( groups[[ig]]$dimid )
				stop("Error, cannot have NAs as dimids!")
				}

			if( verbose )
				print(paste("nc_open: getting dim info for dim number",idim,"in group \"",
					groups[[ig]]$name, "\" dim ID=", dimid2use))

			#-----------------------------------------------------------
			# As a general note, the function ncdim_inq does NOT return
			# a full-fledged "ncdim" object. It returns only a subset of 
			# fields that directly correspond to a low-level, netCDF 
			# dimension in a file.  We now fill in the rest of the fields 
			# to make it into a real ncdim object.
			#-----------------------------------------------------------
			d <- ncdim_inq(groups[[ig]]$id,dimid2use)  # sets dim$name, $len, $unlim

			#---------------------------------------------------------------------
			# Routine ncdim_inq sets only the simple name, not the fully qualified
			# dim name.  Fix that now
			#---------------------------------------------------------------------
			if( groups[[ig]]$name != "" )	# this comparison is FALSE if this is the root group
				d$name <- paste( groups[[ig]]$name, "/", d$name, sep='' )	# example: "model1/run1/longitude".  NO leading slash!

			d$group_index 	<- ig
			d$group_id	<- groups[[ig]]$id
			d$id	   	<- dimid2use				# note: dim$id is the raw C-style integer ID to use WITHIN THE CORRECT GROUP

			#------------------
			# Handle the dimvar
			#------------------
			tt 		<- ncvar_id(groups[[ig]]$id, nc4_basename(d$name))	# note: dimvarid must be used with correct group.  This is -1 if there is no dimvar
			d$dimvarid = ncdf4_make_id( id=tt, group_index=ig, group_id=groups[[ig]]$id, list_index=-1, isdimvar=TRUE )	# NOTE: dimvars are not on the global var list, so list_index == -1

			if( verbose )
				print(paste(".....dim name is",d$name,"  id=", d$id, "  len=",d$len,"     dimvarid=",d$dimvarid$id))

			if( d$dimvarid$id == -1 ) {	# No dimvar for this dim
				d$vals  <- 1:d$len
				d$units <- ""
				d$create_dimvar <- FALSE	# in case this dim is passed to nc_create()
				}
			else {	
				# This dim has a dimvar -- get its properties
				if( verbose )
					print(paste("nc_open: getting dimvar info for dim ",d$name))

				attv <- ncatt_get_inner( d$dimvarid$group_id, d$dimvarid$id, "units" )
				if( attv$hasatt )
					d$units <- attv$value
				else
					d$units <- ""

				attv <- ncatt_get_inner( d$dimvarid$group_id, d$dimvarid$id, "calendar" )
				if( attv$hasatt )
					d$calendar <- attv$value
				#else
				#       since nothing else is defined here, is.null(d$calendar) will return TRUE

				if( d$unlim && (! readunlim)) # Is unlimited, don't read vals, too slow
					d$vals <- rep(NA,d$len)
				else			# Otherwise, read vals
					d$vals <- ncvar_get_inner( d$dimvarid$group_id, d$dimvarid$id, default_missval_ncdf4(), 
							verbose=verbose )

				d$create_dimvar <- TRUE		# in case this dim is passed to nc_create()
				}
			attr(d,"class") <- "ncdim4"	# Is a complete ncdim object now

			global_dim_counter <- global_dim_counter + 1	# global is for ALL groups
			if( verbose )
				{
				print("------------------------------")
				print("Here is new dim:")
				print(paste("Global index=", global_dim_counter, 
					"name=",d$name,"len=",d$len,"unlim=",d$unlim,
					"id=",d$id,"dimvarid=",d$dimvarid$id,"units=",d$units))
				print("------------------------------")
				}

			if( d$unlim && (nc$unlimdimid == -1))	
				nc$unlimdimid = global_dim_counter	# NOTE: technically, this is not an ID, it's an index into the nc$dim[[]] list

			nc$dim[[global_dim_counter]] <- d		# NOTE: nc$dim[[]] list is across ALL dims in file, regardless of group
			dimnames[global_dim_counter] <- d$name		# NOTE: in ncdf4, dim[[]] is indexed by FULLY QUALIFIED dim name
			if( verbose )
				print(paste(".......nc_open: done processing dim ",d$name))
			}
		}

	attr(nc$dim,"names") <- dimnames
	if(verbose) {
		print("nc_open: setting dim$<names> to:")
		print(dimnames)
		}
	
	#-------------------------------------------
	# Get all the vars that this file has.  Note
	# that dimvars are NOT included in the count
	# of vars!!
	#-------------------------------------------
	if( verbose )
		print(paste("nc_open: getting var info.  Number of vars (INCLUDING dimvars)=",tot_nvars_inc_dimvars))
	nc$nvars <- 0			# NOTE this is a GLOBAL list, i.e., it includes vars from ALL groups
	nc$var   <- list()
	varnames <- character()
	have_warned_noncompliant <- FALSE
	for( ig in 1:length(groups)) {

		for( ivar in nc4_loop(1,groups[[ig]]$nvars)) {	

			#-----------------------------------------------------------------
			# Note this is the 'simple' name, NOT the fully qualified var name
			#-----------------------------------------------------------------
			name <- ncvar_name( groups[[ig]]$id, ivar-1 )	# ivar-1 because ncvar_name takes input in C standard, which starts at 0

			if( verbose ) print(paste("Working on group",ig,"(of",length(groups),"), var", ivar, "(of", groups[[ig]]$nvars,"), name=", name))
			if( ncdim_id( groups[[ig]]$id, name ) == -1 ) {	# Only process if NOT a dimvar
				#--------------------------------------
				# No dim with same name as this var, so
				# this var must NOT be a dimvar.
				#--------------------------------------
				if( verbose )
					print(paste("nc_open var loop: will process with group id=", groups[[ig]]$id, " varid=",ivar,"  var name=",name))

				#---------------------------------------------------------------------------------------------------
				# ncvar_inq sets the following:
				#	$id: an object of class 'ncid', with $id, $group_id, $group_index(==-1), and list_index(==-1)
				#	$name
				#	$ndims
				#	$natts
				#	$size
				#	$dimids: raw integer C-style (0-based counting) dimids but in R order
				#	$prec: one of "short', "int", "float", "double", "byte"
				#	$units: units string, or ""
				#	$longname: longname attribute, or ""
				#---------------------------------------------------------------------------------------------------
				v <- ncvar_inq( groups[[ig]]$id, ivar-1 ) # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
				attr(v,"class") <- "ncvar4"
				nc$nvars <- nc$nvars + 1
				v$group_index = ig
				v$id$list_index  = nc$nvars

				#--------------------------------------------------------------------------------
				# Routine ncvar_inq sets the var's name to the SIMPLE (not-fully qualified) name.
				# Fix that now.
				#--------------------------------------------------------------------------------
				if( groups[[ig]]$name != "" )	# this comparison is FALSE if this is the root group
					v$name <- paste( groups[[ig]]$name, "/", v$name, sep='' )	# example: "model1/run1/Temperature".  NO leading slash!

				#---------------------------------------------------
				# Get netcdf-4 specific information for the variable
				#---------------------------------------------------
				if( nc$format == 'NC_FORMAT_NETCDF4' ) {

					#-----------------------
					# Inquire about chunking
					#-----------------------
					chunkrv = ncvar_inq_chunking( groups[[ig]]$id, ivar-1, v$ndims ) # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
					v$chunksizes = chunkrv$chunksizes
					v$storage    = chunkrv$storage		# 1 for NC_CONTIGUOUS and 2 for NC_CHUNKED

					#--------------------------
					# Inquire about compression
					#--------------------------
					comprv = ncvar_inq_deflate( groups[[ig]]$id, ivar-1 ) # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
					v$shuffle = comprv$shuffle
					if( comprv$deflate == 0 )
						v$compression = NA
					else
						v$compression = as.integer( comprv$deflate_level )
					}
				else
					{
					v$chunksizes  = NA
					v$storage     = 1		# 1 for NC_CONTIGUOUS and 2 for NC_CHUNKED
					v$shuffle     = FALSE
					v$compression = NA
					}

				#-----------------------------------------------------------------------------------
				# Get this var's dims.  Special netcdf-4 note: it seems that the way the netcdf-4
				# library works is that variables are numbered starting at 0 for each group, but
				# that dims are numbered consecutively across all groups.  I guess this is mandated
				# by the way that a var's dims can be in either the var's group or ANY PARENT GROUP,
				# which implies that dimids cannot ever be repeated in a file, even if the file has
				# multiple groups.
				#-----------------------------------------------------------------------------------
				v$dims   <- list()
				varunlim <- FALSE
				if( v$ndims > 0 ) {
					for( j in 1:v$ndims ) {
						dimid2find = v$dimids[j]
						matchidx = -1
						for( iidim in 1:nc$ndims ) {
							if( nc$dim[[iidim]]$id == dimid2find ) {
								matchidx = iidim
								break
								}
							}
						if( matchidx == -1 )
							stop(paste("internal error, did not find dim with id=",dimid2find,"in dim list!"))
						v$dim[[j]] = nc$dim[[matchidx]]

						if( v$dim[[j]]$unlim )
							varunlim <- TRUE
						v$varsize <- append(v$varsize, v$dim[[j]]$len)
						}
					}
				v$unlim <- varunlim

				#------------------------------------------------------------
				# Fix for classic files that are unlimited -- they are stored
				# in chunked format, as far as I know
				#------------------------------------------------------------
				if( (nc$format != 'NC_FORMAT_NETCDF4') && varunlim) 
					v$storage = 2		# 1 for NC_CONTIGUOUS and 2 for NC_CHUNKED

				#----------------------------------------
				# Get this var's missing value, or set to
				# a default value if it does not have one
				#----------------------------------------
				found_mv <- FALSE
				v$make_missing_value <- FALSE
				mv <- ncatt_get_inner( groups[[ig]]$id, ivar-1, "missing_value" )  # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
				if( mv$hasatt ) {
					found_mv <- TRUE
					v$missval <- mv$value
					v$make_missing_value <- TRUE
					}
				else
					{
					mv <- ncatt_get_inner( groups[[ig]]$id, ivar-1, "_FillValue" )  # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
					if( mv$hasatt ) {
						found_mv <- TRUE
						v$missval <- mv$value
						v$make_missing_value <- TRUE
						}
					}

				if( ! found_mv ) {
					if( (v$prec=="float") || (v$prec=="double"))
						v$missval <- default_missval_ncdf4()
					else
						v$missval <- NA
					}

				#--------------------------------------------
				# Special check for noncompliant netCDF files
				#--------------------------------------------
				if( (v$prec=="float") || (v$prec=="double")) {
					if( storage.mode(v$missval) == "character" ) {
						v$missval <- as.double( v$missval )
						if(! have_warned_noncompliant ) {
							print(paste("WARNING file",filename,"is not compliant netCDF; variable",name," is numeric but has a character-type missing value! This is an error!  Compensating, but you should fix the file!"))
							have_warned_noncompliant <- TRUE 
							}
						}
					}

				#-------------------------------------------
				# Get add_offset and scale_factor attributes 
				#-------------------------------------------
				ao <- ncatt_get_inner( groups[[ig]]$id, ivar-1, "add_offset" )  # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
				if( ao$hasatt ) {
					v$hasAddOffset <- TRUE
					v$addOffset    <- ao$value
					}
				else
					v$hasAddOffset <- FALSE
				sf <- ncatt_get_inner( groups[[ig]]$id, ivar-1, "scale_factor" )  # ivar-1 because ncvar_inq takes input in C standard, which starts at 0
				if( sf$hasatt ) {
					v$hasScaleFact <- TRUE
					v$scaleFact    <- sf$value
					}
				else
					v$hasScaleFact <- FALSE
				
				nc$var[[nc$nvars]] <- v
				varnames <- append(varnames,v$name)	# NOTE: fully qualified var name, not simple name
				if( verbose ) {
					print("-----------------------")
					print("Here is new var:")
					print(paste("name=",v$name,"  group_id=",v$id$group_id, "  id=",v$id$id,"   ndims=",v$ndims,"   prec=",v$prec))
					print("size=")
					print(v$size)
					print("dimids=")
					print(v$dimids)
					}
				} # end of "if NOT dimvar" clause
			} # end of loop over vars in this group
		} # end of loop over all groups

	attr(nc$var,"names") <- varnames

	if( verbose )
		print(paste("nc_open: leaving for ncid=",nc$id))

	return(nc)
}

#=========================================================================================================
# This is used to change a missing value to a desired value.  It
# is used when a file has not been created on disk yet, or when
# an existing disk file has been opened to be writable.
# 'varid' may be a character string with the var's name,
# or a ncvar class.  
#
ncvar_change_missval <- function( nc, varid, missval ) {

	if( class(nc) != "ncdf4" ) 
		stop("ncvar_change_missval: passed nc NOT of class ncdf4!")

	#------------------------------------------------------
	# Can't do this if the file is on disk and not writable
	#------------------------------------------------------
	if( (nc$filename != "IN-MEMORY") && (! nc$writable))
		stop("ncvar_change_missval: the netcdf file was NOT opened in write mode!")

	idobj <- vobjtovarid4( nc, varid )	# object of type 'ncid', NOT just a simple integer
	if( idobj == -1 ) 
		stop(paste("error: could not find passed variable in the specified netcdf file. Are you sure it's actually in that file?"))
	idx   <- idobj$list_index
	if( idx < 1 )
		stop(paste("Error, did not specify enough information to identify where var is on the global var list. Are you trying to change the missing value for a dimvar? That operation is not supported"))
	nc$var[[idx]]$missval <- missval

	if( nc$filename != "IN-MEMORY" )
		ncatt_put( nc, varid, "missing_value", missval )
}

#=========================================================================================================
# This is the public interface for creating a netCDF file
# on disk.  It takes arguments of class ncvar to put in the file.
# It creates the file on disk and returns an object of class "ncdf4"
# that can be used to access that file.
#
# Example usage, where "temp" and "salin" are objects of class ncvar4:
#
#	nc <- ncdf.create( "test.nc", list(temp,salin))
# or
#	nc <- ncdf.create( "test.nc", salin )
#
nc_create <- function( filename, vars, force_v4=FALSE, verbose=FALSE ) {

	if( verbose ) print(paste('nc_create: entering, package version', nc_version() ))

	if( (! is.character(filename)) || (nchar(filename)<1))
		stop("input filename must be a character string")

	#----------------------------------------------------
	# Have to tell if the input vars is a single var or a 
	# list of vars.   Do it by examining vars$class.  If 
	# vars is a single var, then this will equal "ncvar"; 
	# if vars orig is a list of vars, this will be NULL.
	#----------------------------------------------------
	if( is.character(class(vars)) && (class(vars) == "ncvar4") ) {
		vars <- list(vars)
		if( verbose )
			print("nc_create: input was a single var")
		}
	else if(is.character(class(vars)) && (class(vars) == "list") ) { 
		if( length(vars) < 1 ) 
			stop("Error, at least one ncvar object must be supplied in the vars list")
		if( (! is.character(class(vars[[1]]))) || (class(vars[[1]]) != "ncvar4"))
			stop("Error, second arg must either be a ncvar object (created by a call to ncvar_def()) or a list of ncvar objects")
		#-------------------------------------------------------
		# Make sure ALL elements in the list are of class ncvar4
		#-------------------------------------------------------
		for( ilist in nc4_loop(2,length(vars)) )
			if( (! is.character(class(vars[[ilist]]))) || (class(vars[[ilist]]) != "ncvar4"))
				stop(paste("Error, found an element of the vars list that is NOT an object created by a call to ncvar_def...element #",ilist,sep=''))
		if( verbose )
			print("nc_create: input was a list of vars")
		}
	else
		stop("Error, second arg must either be a ncvar object (created by a call to ncvar_def()) or a list of ncvar objects")

	#------------------------------------------------------------
	# Figure out the list of all groups that we will be creating.
	#------------------------------------------------------------
	if( verbose ) print('nc_create: parsing group structure')
	group <- nc_parse_group_structure( vars, verbose=verbose )
	if( length(group) > 1 ) {
		if( verbose ) print("Forcing netcdf version 4 format file since there is more than 1 group")
		force_v4 <- TRUE
		}

	#-----------------------------------------------------------------------
	# If any variables use compression or chunking, we must create a V4 file
	#-----------------------------------------------------------------------
	if( verbose ) print('nc_create: checking to see if we MUST produce a netcdf-version 4 file')
	for( ivar in 1:length(vars)) {
		use_shuffle     = vars[[ivar]]$shuffle
		use_compression = (! is.na(vars[[ivar]]$compression))
		use_chunking    = ((length(vars[[ivar]]$chunksizes)>1) || (! is.na(vars[[ivar]]$chunksizes)))
		if( verbose ) print(paste("var: >", vars[[ivar]]$name,"<  Use shuffle: >",use_shuffle,
				"<  use_compression: >",use_compression,
				"<  use_chunking: >",use_chunking,"<", sep='' ))
		if( use_shuffle || use_compression || use_chunking ) {
			if( verbose ) print("Forcing netcdf version 4 format file since a var using compression or chunking")
			force_v4 = TRUE
			}
		else
			{
			if( verbose ) print("Not forcing netcdf version 4 format file since no var is using compression or chunking")
			}
		}

	#---------------------------------------------------------------------
	# If there are multiple unlimited dimensions, we must create a V4 file
	#---------------------------------------------------------------------
	if( verbose ) print('nc_create: checking to see if there are multiple unlimited dims (which would force V4 file)')
	unlim_dimname = ''
	multi_unlim_dims = FALSE
	for( ivar in 1:length(vars)) {
		ndims = vars[[ivar]]$ndims
		if( ndims > 0 ) {
			for( idim in 1:ndims ) {
				if( vars[[ivar]]$dim[[idim]]$unlim ) {
					if( nchar(unlim_dimname) == 0 ) 
						unlim_dimname = vars[[ivar]]$dim[[idim]]$name
					else
						if( vars[[ivar]]$dim[[idim]]$name != unlim_dimname ) {
							if( verbose ) print(paste("Forcing netcdf version 4 format file",
								"since there is more than 1 unlimited dim"))
							force_v4 = TRUE
							multi_unlim_dims = TRUE
							break
							}
					}
				}
			}
		}
	if( verbose ) {
		if( multi_unlim_dims )
			print('nc_create: Yes, there ARE multiple unlimited dims in this file, forcing V4')
		else
			print('nc_create: No, there are not any multiple unlimited dims in this file')
		}

	nc <- list()

	#-----------------------------------------------------------
	# These values MUST MATCH the values in the source code file 
	# ncdf.c, routine R_nc4_create!
	#-----------------------------------------------------------
	flag_NC_NOCLOBBER 	<- 1
	flag_NC_SHARE     	<- 2
	flag_NC_64BIT_OFFSET	<- 4
	flag_NC_NETCDF4		<- 8

	#----------------
	# Create the file
	#----------------
	nc$cmode    <- 0
	if( force_v4 )
		nc$cmode <- nc$cmode + flag_NC_NETCDF4
	nc$error    <- -1
	nc$id       <- -1
	if( verbose )
		print(paste("Calling R_nc4_create for file ",filename))
	nc<-.C("R_nc4_create",
		filename,
		as.integer(nc$cmode),
		id=as.integer(nc$id),
		error=as.integer(nc$error),
		PACKAGE="pbdNCDF4")
	if( nc$error != 0 )
		stop("Error in nc_create!")
	if( verbose )
		print(paste("back from R_nc4_create for file ",filename))
	nc$nvars  <- 0
	attr(nc,"class")  <- "ncdf4"
	nc$filename <- filename
	nc$writable <- TRUE

	nc$ndims  <- 0
	nc$dim    <- list()
	nc$var    <- list()

	#--------------------------
	# Create groups in the file
	#--------------------------
	nc$group 	<- group
	nc$ngroups    	<- length( group )	# Note: will always be at least 1, since root group is in there
	nc$group[[1]]$id <- nc$id		# id of group 1, which is the root group, is always exactly equal to ncid
	nc$fqgn2Rindex 	<- list()		# given fqgn, gives R index into nc$group list
	nc$fqgn2Rindex[["/"]] <- 1		# the root group is always the first entry on the group list
	if( nc$ngroups > 1 ) {
		if( verbose )
			print(paste("nc_create: about to create the",nc$ngroups-1,"non-root groups"))
		for( ig in 2:nc$ngroups ) {	# Note: skip root group
			group_name <- nc$group[[ig]]$name
			fqgn       <- nc$group[[ig]]$fqgn		# fully qualified group name
			fqpn       <- nc$group[[ig]]$fqpn		# fully qualified parent group name
			if( fqpn == "" ) fqpn <- '/'
			pidx       <- nc$fqgn2Rindex[[ fqpn ]]		# index of parent group in group list
			if( verbose )
				print(paste("Defining group", group_name, "(which has parent", fqpn,")" ))
			if( (pidx<1) || (pidx>nc$ngroups) || is.null(pidx)) {
				print(paste("error, did not find index for parent group", fqpn, "in fqgn2Rindex!"))
				print('Here are the entries in fqgn2Rindex:')
				nn <- length(nc$fqgn2Rindex)
				for( ii in 1:nn )
					print(paste(names(nc$fqgn2Rindex)[ii], '-->', nc$fqgn2Rindex[[ii]] ))
				stop("internal error")
				}
			parentid  <- nc$group[[pidx]]$id
			if( parentid < 1 ) 
				stop(paste("error, apparently the parent group", fqpn, "has not been defined yet, but i need to reference it!"))

			#-------------------------------------------
			# This call actually makes the group on disk
			#-------------------------------------------
			gid <- nc_make_group_inner( parentid, group_name )
			if( verbose )
				print(paste("back from nc_make_group_inner call for file ",filename,"; just made group ", group_name," with new group id:", gid))

			nc$group[[ig]]$id <- gid
			nc$fqgn2Rindex[[ fqgn ]] <- ig
			}
		}

	max_nc_dims <- 20

	#---------------------------------------------------------------
	# Add the vars to the file.  NOTE that this also adds the unique
	# dims (and hence, dimvars) to the file as a side effect!
	# Note also that the returned 'nc' value is updated each time
	# this subroutine is called.
	#---------------------------------------------------------------
	if( verbose ) print("nc_create: about to create the vars")
	for(ivar in 1:length(vars)) 
		nc <- ncvar_add( nc, vars[[ivar]], verbose=verbose, indefine=TRUE )

	#-----------------------------------------------------------
	# Set the names attribute on the $var and $dim lists so that
	# we can access them by FULLY QUALIFIED name instead of only 
	# by position
	#-----------------------------------------------------------
	if( verbose ) print("nc_create: setting names attribute on $var and $dim lists")
	varnames <- array('',nc$nvars)
	for( ivar in nc4_loop(1,nc$nvars))
		varnames[ivar] <- nc$var[[ivar]]$name
	attr(nc$var,"names") <- varnames
	dimnames <- array('',nc$ndims)
	for( idim in nc4_loop(1,nc$ndims) )
		dimnames[idim] <- nc$dim[[idim]]$name
	attr(nc$dim,"names") <- dimnames

	#-----------------
	# Exit define mode
	#-----------------
	if( verbose ) print("nc_create: exiting define mode")
	nc_enddef( nc )

	#-----------------------------------------------------------------------
	# Have to set the format string of the ncdf4 object.  We could calculate
	# this from first principles by analyzing the kind of file we just 
	# created.  However I'm thinking it's easier and less error prone simply
	# to read back in the format from the file that's now on disk.  
	# Possible (string) values are:
	# 'NC_FORMAT_CLASSIC', 'NC_FORMAT_64BIT', 'NC_FORMAT_NETCDF4',
	# and 'NC_FORMAT_NETCDF4_CLASSIC'
	#-----------------------------------------------------------------------
	nc$format = ncdf4_format( nc$id )
	if( verbose )
		print(paste("file", filename, "is format", nc$format ))

	return(nc)
}

#===============================================================
# This is a SPECIAL PURPOSE function ONLY to be used when adding
# an already defined variable (accomplished via "ncvar_def")
# to an ALREADY EXISTING netcdf file.  Normally, when making
# a new netcdf file from scratch, a list of vars to be created
# would be passed to "nc_create"; this is the preferred method
# of putting vars in a file.  However, sometimes it's necessary
# to add a var to an already-existing file; that's what this
# routine is for.
#
# Inputs:
#	nc : object of class 'ncdf4'
#	v  : object of class 'ncvar4'
#
ncvar_add <- function( nc, v, verbose=FALSE, indefine=FALSE ) {
	
	if( verbose )
		print(paste("ncvar_add: entering with indefine=",indefine))

	if( class(nc) != "ncdf4" ) 
		stop("ncvar_add: passed nc NOT of class ncdf4!")
	if( verbose )
		print(paste("ncvar_add: ncid of file to add to=",nc$id,
			"   filename=",nc$filename,"    writable=",nc$writable))

	if( class(v) != "ncvar4" ) 
		stop("var.add.ncdf: passed var NOT of class ncvar4! The second arg to var.add.ncdf must be the return value from a call to ncvar_def")
	if( verbose )
		print(paste("ncvar_add: varname to add=",v$name))

	if( ! indefine ) {
		if( verbose )
			print(paste("ncvar_add: about to redef ncid=",nc$id))
		nc_redef( nc )	# Go back into define mode
		}

	#-----------------------------------------------------
	# Create the dims for this var.  Harder than it sounds 
	# because we must take care not to repeat making a dim 
	# that occurs in more than one variable.
	#---------------------------------------------------
	nd <- v$ndims
	dimvarids <- array(0,nd)
	if( verbose )
		print(paste("ncvar_add: creating",nd,"dims for var",v$name))
	for( idim in nc4_loop(1,nd) ) {
		d <- v$dim[[idim]]
		if( verbose )
			print(paste("ncvar_add: working on dim >",d$name,"< (number",idim,") for var",v$name))

		#-----------------------------------------------
		# See if we've already made a dim with this name
		# in this group
		#-----------------------------------------------
		place <- -1
		for( ii in nc4_loop(1,length(nc$dim))) {
			if( nc$dim[[ii]]$name == d$name ) {
				#---------------------------------------------------------------
				# Check to make sure this is REALLY the same dim, even though we
				# know it has the same name as an existing dim!
				#---------------------------------------------------------------
				if( ! ncdim_same( nc$dim[[ii]], d )) {
					stop(paste("Error, when trying to add variable named",
						v$name, "to file",nc$filename,"I found this variable has a dim named",d$name,
						"However, the file ALREADY has a dim named",nc$dim[[ii]]$name,
						"with different characteristics than the new dim with the same name!",
						"This is not allowed."))
					}
				if( verbose )
					print(paste("ncvar_add: dim",d$name, "has been seen before"))
				place <- ii
				break
				}
			}
		if( place == -1 ) {
			#--------------------------------------------
			# This dim has not been seen before -- create
			#--------------------------------------------
			if( verbose )
				print(paste("ncvar_add: creating dim",d$name, "(which has not been seen before)"))

			#--------------------------------------------------------------------------------
			# Following call returns dim ID, divar ID, group ID used to access the first two,
			# and group INDEX used to access the first two
			#--------------------------------------------------------------------------------
			ids         	<- ncdim_create(nc,d,verbose)	# *** NOTE: makes the dimvar, too! ***
			dimid       	<- ids[1]
			dimvarid    	<- ids[2]
			group_id 	<- ids[3]
			group_index 	<- ids[4]

			newel       	<- list()
			attr(newel,"class") <- "ncdim4"
			newel$name     	<- d$name
			newel$units    	<- d$units
			newel$vals     	<- d$vals
			newel$len     	<- d$len
			newel$unlim    	<- d$unlim

			nc$ndims    	<- nc$ndims + 1
			newel$id       	<- dimid		# remember, dimids are just simple integers, varids are ncid objects
			newel$dimvarid 	<- ncdf4_make_id( id=dimvarid, group_index=group_index, 
					group_id=group_id, list_index=nc$ndims, isdimvar=TRUE )
			nc$dim[[nc$ndims]] <- newel
			}
		else
			dimid <- nc$dim[[place]]$id	# simple C-style 0-counting based integer

		dimvarids[idim] <- dimid	# OK that this is the simple integer, not a ncid object, 
						# since by v4 spec dims must be visible to vars.  
		}

	#----------------------------------------------------
	# Reverse the dimvarids, because R uses Fortran-style
	# ordering and we are using the C netCDF interface
	#----------------------------------------------------
	dimids <- dimvarids
	if( nd > 0 ) 
		dimids <- dimids[length(dimids):1]
	newvar       <- list()
	newvar$id    <- -1
	newvar$error <- -1

	#-----------------------------------------------------------
	# Select the routine we will be using to create the variable
	# based on its precision
	#-----------------------------------------------------------
	if( verbose )
		print(paste("ncvar_add: creating",v$prec,"precision var",v$name))
	if( (v$prec == "integer") || (v$prec == "int") )
		funcname <- "R_nc4_def_var_int"
	else if( v$prec == "short" )
		funcname <- "R_nc4_def_var_short"
	else if( (v$prec == "float"))
		funcname <- "R_nc4_def_var_float"
	else if( v$prec == "double" )
		funcname <- "R_nc4_def_var_double"
	else if( v$prec == "char" )
		funcname <- "R_nc4_def_var_char"
	else if( v$prec == "byte" )
		funcname <- "R_nc4_def_var_byte"
	else
		stop(paste("internal error in nc_create: var has unknown precision:",v$prec,". Known vals: short float double integer char byte"))

	#-----------------------------------------------------------------------
	# Figure out the ncid to use.  If this var is in the root group, it will
	# just be the ncid.  Otherwise, it will be the group id
	#-----------------------------------------------------------------------
	if( verbose ) print('ncvar_add: figuring out ncid to use')
	if( nslashes_ncdf4( v$name ) == 0 ) 
		gidx <- 1
	else
		{
		vars_fqgn <- nc4_basename( v$name, dir=TRUE )	# this is the var's fully qualified GROUP name
		gidx      <- nc$fqgn2Rindex[[ vars_fqgn ]]
		if( is.null(gidx))
			stop(paste("internal error: did not find fully qualified group name", vars_fqgn," in list of groups for file", nc$filename))
		}
	ncid2use <- nc$group[[gidx]]$id
	name2use <- nc4_basename( v$name )
	if( verbose ) print(paste('ncvar_add: ncid2use=', ncid2use, 'name2use=', name2use ))

	#---------------------------------
	# Now actually create the variable
	#---------------------------------
	### WCC: R CMD check warning.
	# newvar<-.C(funcname,
	# 	as.integer(ncid2use),
	# 	as.character(name2use),
	# 	as.integer(v$ndims),
	# 	as.integer(dimids),	
	# 	id=as.integer(newvar$id),
	# 	error=as.integer(newvar$error),
	# 	PACKAGE="pbdNCDF4")
	### WCC: avoid R CMD check warning.
	if( (v$prec == "integer") || (v$prec == "int") ){
		newvar<-.C("R_nc4_def_var_int",
			as.integer(ncid2use),
			as.character(name2use),
			as.integer(v$ndims),
			as.integer(dimids),	
			id=as.integer(newvar$id),
			error=as.integer(newvar$error),
			PACKAGE="pbdNCDF4")
	} else if( v$prec == "short" ){
		newvar<-.C("R_nc4_def_var_short",
			as.integer(ncid2use),
			as.character(name2use),
			as.integer(v$ndims),
			as.integer(dimids),	
			id=as.integer(newvar$id),
			error=as.integer(newvar$error),
			PACKAGE="pbdNCDF4")
	} else if( v$prec == "float" ){
		newvar<-.C("R_nc4_def_var_float",
			as.integer(ncid2use),
			as.character(name2use),
			as.integer(v$ndims),
			as.integer(dimids),	
			id=as.integer(newvar$id),
			error=as.integer(newvar$error),
			PACKAGE="pbdNCDF4")
	} else if( v$prec == "double" ){
		newvar<-.C("R_nc4_def_var_double",
			as.integer(ncid2use),
			as.character(name2use),
			as.integer(v$ndims),
			as.integer(dimids),	
			id=as.integer(newvar$id),
			error=as.integer(newvar$error),
			PACKAGE="pbdNCDF4")
	} else if( v$prec == "char" ){
		newvar<-.C("R_nc4_def_var_char",
			as.integer(ncid2use),
			as.character(name2use),
			as.integer(v$ndims),
			as.integer(dimids),	
			id=as.integer(newvar$id),
			error=as.integer(newvar$error),
			PACKAGE="pbdNCDF4")
	} else if( v$prec == "byte" ){
		newvar<-.C("R_nc4_def_var_byte",
			as.integer(ncid2use),
			as.character(name2use),
			as.integer(v$ndims),
			as.integer(dimids),	
			id=as.integer(newvar$id),
			error=as.integer(newvar$error),
			PACKAGE="pbdNCDF4")
	} else{
		stop(paste("internal error in nc_create: var has unknown precision:",v$prec,". Known vals: short float double integer char byte"))
	}
	if( verbose )
		print(paste("nc_create: C call returned value",newvar$error))
	if( newvar$error != 0 ) {
		print('----------------------')
		print(paste('Var: ', v$name))
		print(paste('Ndims: ', v$ndims))
		print('Dimids: ')
		print(dimids)
		stop(paste("Error in ncvar_add, defining var", v$name))
		}

	nc$nvars   <- nc$nvars + 1
	v$id       <- ncdf4_make_id( newvar$id, gidx, ncid2use, nc$nvars )
	nc$var[[nc$nvars]] <- v

	#----------------------------------------
	# Set compression parameters if requested
	#----------------------------------------
	if( v$shuffle || (! is.na(v$compression))) {
		if( verbose ) print(paste("nc_var_add: setting compression params"))
		if( v$shuffle )
			shuffle_param = 1
		else
			shuffle_param = 0
		if( is.na(v$compression) ) {
			d_param = 0
			d_level = 0
			}
		else
			{
			d_param = 1
			d_level = v$compression
			}
		ncvar_def_deflate( ncid2use, newvar$id, shuffle_param, d_param, d_level )
		}

	#-------------------------------------
	# Set chunking parameters if requested
	#-------------------------------------
	if( (length(v$chunksizes)>1) || (! is.na(v$chunksizes)) ) {
		if( verbose ) print(paste("nc_var_add: setting chunking params"))
		chunksizes = as.integer(v$chunksizes)
		#--------------------------------------------------------------
		# Make sure no chunksize is larger than a fixed (non-unlim) dim
		#--------------------------------------------------------------
		if( length(chunksizes) != v$ndims )
			stop(paste("Error, for var",v$name,"ndims=",v$ndims,"but length of chunksizes array=",
				length(chunksizes),".  They must be the same!"))
		for( ii in 1:v$ndims ) {
			if( (! v$dim[[ii]]$unlim) && (chunksizes[ii] > v$dim[[ii]]$len))
				stop(paste("Error in supplied chunksizes: dim number",ii,
					", named", v$dim[[ii]]$name,", is length",
					v$dim[[ii]]$len,"but chunksizes for this dim is ",
					chunksizes[ii],". Chunksizes must be <= dim length!"))
			if( chunksizes[ii] < 1 )
				stop(paste("Error in supplied chunksizes: dim number",ii,
					", named", v$dim[[ii]]$name,", is length",
					v$dim[[ii]]$len,"but chunksizes must be >= 1"))
			}
		storage = 2	# set flag to enable chunking
		ncvar_def_chunking( ncid2use, newvar$id, storage, chunksizes )
		}

	#------------------------------------------------------
	# Add the attributes -- units, missing_value, long_name
	#------------------------------------------------------
	if( verbose ) print(paste("nc_var_add: adding attributes"))
	if( (! is.null( v$units )) && (! is.na(v$units)) && (nchar(v$units)>0)) 
		ncatt_put_inner( ncid2use, newvar$id, "units", v$units, definemode=TRUE )

	#------------------------------------
	# Create a missing value if requested
	#------------------------------------
	if( verbose ) print(paste("nc_var_add: creating missing value if requested"))
	if( v$make_missing_value ) {
		if( is.null( v$missval ) && ((v$prec=="float") || (v$prec=="double"))) {
			ncatt_put_inner( ncid2use, newvar$id, "_FillValue", default_missval_ncdf4(), definemode=TRUE )
			}
		else
			{
			if( ! is.null(v$missval) ) {
				ncatt_put_inner( ncid2use, newvar$id, "_FillValue", v$missval, definemode=TRUE, verbose=verbose, prec=v$prec )
				}
			}
		}

	#-------------------------
	# Make long name attribute
	#-------------------------
	if( verbose ) print(paste("nc_var_add: creating long name attribute"))
	if( (v$longname != v$name) && (nchar(v$longname)>0))
		ncatt_put_inner( ncid2use, newvar$id, "long_name", v$longname, definemode=TRUE )

	if( ! indefine )
		nc_enddef( nc )	# Exit define mode

	if( verbose ) print(paste("nc_var_add: returning"))
	return(nc)
}

#===============================================================
# This has two modes, depending on the value of 'attname' that
# is passed in.
#
# If an attname is specified, then 
# this returns a list; first element of list (named "hasatt")
# is TRUE if the variable had an attribute with name "attname", and
# is FALSE otherwise.  Second element of the list (named "value")
# holds the value IFF the variable had an attribute of that name.
# "varid" can be an integer (R-type varid with counting starting
# at 1), object of class ncvar4, or character string with a var's name.
# If the attribute type is short or integer, an integer value is
# returned.  If the attribute type is float or double, a double
# value is returned.  If the attribute type is text, a character
# value is returned.
#
# If no attname is specified, then this returns a list, with the
# name of the element set to the attribute name, and the value of
# the element set to the attribute's value.  If the varid is 0,
# then global attributes are returned.  If the indicated variable
# (or file) has no attributes, then a list with 0 elements is
# returned.
#
ncatt_get <- function( nc, varid, attname=NA, verbose=FALSE ) {

	if( verbose ) print('ncatt_get: entering')

	if( class(nc) != "ncdf4" ) 
		stop("Error, first passed argument must be an object of class ncdf4")

	#----------------------------------------------------------------------------
	# Atts have a special case where an integer 0 means to access the global atts
	#----------------------------------------------------------------------------
	if( is.numeric(varid) && (varid == 0)) {
		is_global = TRUE
		if( verbose ) print('ncatt_get: is a global att')
		}
	else
		{
		if( (class(varid) != "ncvar4") && (!is.character(varid)))
			stop(paste("second arg (varid) must be one of: 0, an object of class ncvar4, or the character string name of a variable"))
		if( verbose ) print('ncatt_get: is NOT a global att')
		is_global = FALSE
		}

	if( is_global ) {
		if( verbose ) print('ncatt_get: calling ncatt_get_inner for a global att')
		return( ncatt_get_inner( nc$id, -1, attname=attname, verbose=verbose ))	# NOTE how we convert from varid=0 to varid=-1 here to match C API standard
		}
	else
		{
		if( (class(varid) != "ncvar4") && ( storage.mode(varid) != "character" )) 
			stop("ncvar_change_missval: error, passed varid must be either 0 (for global attributes), the name of the variable to operate on, or an object of class ncvar4")
		if( verbose ) print('ncatt_get: getting object id')
		idobj <- vobjtovarid4( nc, varid, allowdimvar=TRUE, verbose=verbose )	# an object of class 'ncid4'
		if( idobj$id == -1 )
			return( list() )
		if( verbose ) print('ncatt_get: calling ncatt_get_inner for a non-global att')
		return( ncatt_get_inner( idobj$group_id, idobj$id, attname=attname, verbose=verbose ))
		}
}

#=================================================================================================
# Put an attribute into a netCDF file.  This puts the file into
# define mode, then takes it back out of define mode when done (unless
# definemode=TRUE, in which case the file is ASSUMED to be in define
# mode already, and is left in define mode as well).
#
# NOTE this does not work with IN-MEMORY ncdf files.
# (Note: you could extend this to work with in-memory files by
# having a list of variable (or file) attributes in the appropriate
# R classes, but I haven't bothered to do this since it does not
# come up very much for me.)
#
# "varid" can be an object of class ncvar4 or character string with a var's name.
# As a special case, if varid==0, then a global attribute is set instead of a variable's
# attribute.
#
# "nc" can be a ncdf4 class object
#
# Precision: ordinarily the precision (type) of the attribute will always
# follow the precision (type) of the var that this is an attribute of.
# However, you can explicitly override this by setting "prec" to the
# desired precision, which can be one of: short float double text( or character) integer.
# In the event of a global attribute, which of course has no associated
# variable, the storage mode of the passed attval will be used to
# determine the precision of the attribute, UNLESS "prec" is set.
# If "prec" is set, it always determines the created attribute type.
#
ncatt_put <- function( nc, varid, attname, attval, prec=NA, 
				verbose=FALSE, definemode=FALSE ) {

	if( verbose ) print('ncatt_put: entering' )

	if( class(nc) != "ncdf4" ) 
		stop("Error, first passed argument must be an object of class ncdf4")

	#-------------------------------------------------------
	# Can't do this if the file is in memory or not writable
	#-------------------------------------------------------
	if( (nc$filename == "IN-MEMORY") || (! nc$writable))
		stop("ncatt_put: the netcdf file has not been written to disk yet, or was not opened in write mode!")

	#----------------------------------------------------------------------------
	# Atts have a special case where an integer 0 means to access the global atts
	#----------------------------------------------------------------------------
	if( verbose ) print('ncatt_put: checking for a global att' )
	if( is.numeric(varid) && (varid == 0)) 
		is_global = TRUE
	else
		{
		if( (class(varid) != "ncvar4") && (!is.character(varid)))
			stop(paste("second arg (varid) must be one of: 0, an object of class ncvar4, or the character string name of a variable"))
		is_global = FALSE
		}

	if( is_global ) {
		if( verbose ) print('ncatt_put: IS a global att' )
		ncatt_put_inner( nc$id, -1, attname, attval, prec=prec, verbose=verbose, definemode=definemode )
		}
	else
		{
		if( verbose ) print('ncatt_put: is NOT a global att' )

		idobj <- vobjtovarid4( nc, varid, allowdimvar=TRUE )	# an object of type "ncid", NOT just a simple integer
		if( verbose ) 
			print(paste("Making attribute",attname,"with value",attval,"for ncid=",idobj$group_id,
				"and varid=",idobj$id))

		#----------------------------------------------------------------------
		# It is possible to identify a dimvar but not actually have that dimvar
		# present.  This is signaled by $isdimvar && $id == -1.
		#----------------------------------------------------------------------
		if( idobj$isdimvar && (idobj$id == -1))
			stop(paste("dimension", varid$name, "in file", nc$filename, "does NOT have a dimvar, so you cannot call ncatt_put with the name of that dimension to try to put an attribute on the dimvar"))

		if( verbose ) print('ncatt_put: calling ncatt_put_inner')
		ncatt_put_inner( idobj$group_id, idobj$id, attname, attval, prec=prec, verbose=verbose, definemode=definemode )
		if( verbose ) print('ncatt_put: back from ncatt_put_inner')
		}

	if( verbose ) print('ncatt_put: exiting')
}

#===============================================================================================
# Writes a vector of values to a netCDF file.  'start' and 'count'
# are given in R convention, i.e., starting at 1, and
# in XYZT order.  If they are omitted, the full array is written.
# 'varid' can be the variable's name or an object of class ncvar4.
# If varid is NA, then the "only" var in 
# the file is assumed to be the selected one.  Values that are NA's 
# in the input data are converted to that variable's 'missing_value' 
# attribute before being written out to the file.
#
# Changs for netcdf-4: if varid is a ncvar object, then all is well.
# Otherwise, if varid is a character string, it must be the fully
# qualified var name.  (Note that it could also be a DIMVAR name.)
#
ncvar_put <- function( nc, varid=NA, vals=NULL, start=NA, count=NA, verbose=FALSE ) {

	if( class(nc) != 'ncdf4' )
		stop(paste("Error: first argument to ncvar_put must be an object of type ncdf,",
			"as returned by a call to nc_open(...,write=TRUE) or nc_create"))

	if( (mode(varid) != 'character') && (class(varid) != 'ncvar4') && (class(varid) != 'ncdim4') && (! is.na(varid)))
		stop(paste("Error: second argument to ncvar_put must be either an object of type ncvar,",
			"as returned by a call to ncvar_def, or the character-string name of a variable",
			"in the file.  If there are multiple vars in the file with the same name (but",
			"in different groups), then the fully qualified var name must be given, for",
			"example, model1/run5/Temperature"))

	#-----------------------------------------------------------------------------------
	# Exactly why we do the following is obscure.  Note that DIMVARS are not kept on 
	# the 'variable' list for the file.  So if we explicitly make a dimvar, then pass
	# that ncvar object to this routine, it will ordinarily fail because no 'var' 
	# matching that dimvar will be found.  This however works if we pass the NAME of
	# the dimvar, because that is matched on the dimvar list as well as the var list.
	# To avoid this error, we force all matches to be by name, rather than by var object
	#-----------------------------------------------------------------------------------
	if( (class(varid) == 'ncvar4') || (class(varid) == 'ncdim4')) {
		varid = varid$name
		if( verbose ) print(paste("ncvar_put: converting passed ncvar4/ncdim4 object to the name:", varid))
		}

	if( is.null(vals))
		stop("requires a vals argument to be set")

	if( verbose ) {
		if( mode(varid) == 'character')
			vname <- varid
		else
			vname <- varid$name
		print(paste("ncvar_put: entering, filename=", nc$filename, ' varname=', vname ))
		}

	#-----------------------------------------------------------------
	# Identify exactly what var (or dimvar) we will be putting data to
	#-----------------------------------------------------------------
	idobj = vobjtovarid4( nc, varid, allowdimvar=TRUE, verbose=verbose )	# NOTE: not a simple integer, but a ncid4 class object with $id, $group_index, $group_id, $list_index
	ncid2use   = idobj$group_id
	varid2use  = idobj$id
	varidx2use = idobj$list_index	# this is the index into the nc$vars[[]] list that indictes this variable
	isdimvar   = idobj$isdimvar

	if( verbose ) 
		print(paste('ncvar_put: writing to var (or dimvar) with id=',idobj$id, ' group_id=', idobj$group_id ))

	#-----------------------------
	# Check inputs for correctness
	#-----------------------------
	sm <- storage.mode(start)
	if( (sm != "double") && (sm != "integer") && (sm != "logical"))
		stop(paste("passed a start argument of storage mode",sm,"; can only handle double or integer"))
	sm <- storage.mode(count)
	if( (sm != "double") && (sm != "integer") && (sm != "logical"))
		stop(paste("passed a 'count' argument with storage mode '",sm,"'; can only handle double or integer", sep=''))

	#--------------------
	# Prevent dumb errors
	#--------------------
	if( ! nc$writable ) 
		stop(paste("trying to write to file",nc$filename,"but it was not opened with write=TRUE"))

	varsize <- ncvar_size ( ncid2use, varid2use )
	ndims   <- ncvar_ndims( ncid2use, varid2use )
	is_scalar = (varsize == 1) && (ndims == 0)
	if( verbose ) {
		print(paste("ncvar_put: varsize="))
		print(varsize)
		print(paste("ncvar_put: ndims=", ndims))
		print(paste("ncvar_put: is_scalar=", is_scalar ))
		}

	#--------------------------------------------------------
	# Fix up start and count to use (in R convention for now)
	#--------------------------------------------------------
	if( (length(start)==1) && is.na(start) ) {
		if( is_scalar )
			start <- 1
		else
			start <- rep(1,ndims)	# Note: use R convention for now
		}
	else
		{
		if( length(start) != ndims ) 
			stop(paste("'start' should specify",ndims,
				"dims but actually specifies",length(start)))
		}
	if( verbose ) {
		print("ncvar_put: using start=")
		print(start)
		}
	if( (length(count)==1) && is.na(count)) {
		count <- varsize - start + 1	
		}
	else
		{
		if( length(count) != ndims ) 
			stop(paste("'count' should specify",ndims,
				"dims but actually specifies",length(count)))
		count <- ifelse( (count == -1), varsize-start+1, count)
		}
	if( verbose ) {
		print("ncvar_put: using count=")
		print(count)
		}

	#------------------------------
	# Switch from R to C convention
	#------------------------------
	c.start <- start[ ndims:1 ] - 1
	c.count <- count[ ndims:1 ]

	#--------------------------------------------
	# Change NA's to the variable's missing value
	#--------------------------------------------
	if( verbose )
		print("about to change NAs to variables missing value")
	if( isdimvar )
		mv <- default_missval_ncdf4()
	else
		mv <- nc$var[[ varidx2use ]]$missval 
	vals <- ifelse( is.na(vals), mv, vals)

	#---------------------------------
	# Get the correct type of variable
	#---------------------------------
	precint <- ncvar_type( ncid2use, varid2use ) # 1=short, 2=int, 3=float, 4=double, 5=char, 6=byte, 7=ubyte, 8=ushort, 9=uint, 10=int64, 11=uint64, 12=string
	if( verbose )
		print(paste("Putting var of type",precint," (1=short, 2=int, 3=float, 4=double, 5=char, 6=byte, 7=ubyte, 8=ushort, 9=uint, 10=int64, 11=uint64, 12=string)"))

	#----------------------------------------------------------
	# Sanity check to make sure we have at least as many values 
	# in the data array as we are writing.  Chars are a special
	# case because typically they are defined with an extra
	# "nchar" dim that is not included in the passed array.
	#----------------------------------------------------------
	n2write <- prod(count)
	if( (precint != 5) && (length(vals) != n2write)) {
		if( length(vals) > n2write ) 
			print(paste("ncvar_put: warning: you asked to write",n2write,
				"values, but the passed data array has",length(vals),
				"entries!"))
		else
			stop(paste("ncvar_put: error: you asked to write",n2write,
				"values, but the passed data array only has",length(vals),
				"entries!"))
		}

	rv <- list()
	rv$error <- -1

	if( verbose ) {
		print("ncvar_put: calling C routines with C-style count=")
		print(c.count)
		print("and C-style start=")
		print(c.start)
		}
	if( (precint == 1) || (precint == 2) || (precint == 6) || (precint == 7) || (precint == 8) || (precint == 9)) {
		#--------------------------------------
		# Short, Int, Byte, UByte, UShort, UInt 
		#--------------------------------------
		rv <- .C("R_nc4_put_vara_int", 
			as.integer(ncid2use),
			as.integer(varid2use),	
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.integer(vals),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		if( rv$error != 0 ) 
			stop("C function R_nc4_put_var_int returned error")
		if( verbose )
			print(paste("C function R_nc4_put_var_int returned", rv$error))
		}

	else if( (precint == 3) || (precint == 4) || (precint == 10) || (precint == 11)) {
		#-----------------------------------------------
		# Float, double, 8-byte int, unsigned 8-byte int
		#-----------------------------------------------
		if( (precint == 10) || (precint == 11)) {
			print(paste(">>>> WARNING <<<< You are attempting to write data to a 8-byte integer,"))
			print(paste("but R does not have an 8-byte integer type.  This is a bad idea! I will"))
			print(paste("TRY to write this by converting from double precision floating point, but"))
			print(paste("this could lose precision in your data!"))
			}
		rv <- .C("R_nc4_put_vara_double", 
			as.integer(ncid2use),
			as.integer(varid2use),	
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.double(vals),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4",
			NAOK=TRUE )
		if( rv$error != 0 ) 
			stop("C function R_nc4_put_var_double returned error")
		if( verbose )
			print(paste("C function R_nc4_put_var_double returned", rv$error))
		}

	else if( precint == 5 ) {
		#----------
		# Character
		#----------
		rv <- .C("R_nc4_put_vara_text", 
			as.integer(ncid2use),
			as.integer(varid2use),	
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.character(vals),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		if( rv$error != 0 ) 
			stop("C function R_nc4_put_var_double returned error")
		if( verbose )
			print(paste("C function R_nc4_put_var_text returned", rv$error))
		}

	else
		stop(paste("Internal error in ncvar_put: unhandled variable type=",precint,". Types I know: 1=short 2=int 3=float 4=double 5=char"))
}

#===============================================================
# Returns data values from a netCDF file.  'start' and 'count'
# are given in R convention, i.e., starting at 1, and
# in XYZT order.  If they are omitted, the full array is read in.
# 'varid' can be the variable's name, an object of class ncvar4,
# or the integer varid.  'varid' can also be omitted entirely,
# in which case the only variable in the file is identified and
# that one read in (if no such variable can be identified, an
# error is generated).  The 'count' array can have -1's, which 
# indicate that all the values in that dim are to be read (subject
# to the start array).  Missing values in the source file (i.e.,
# values that match that variable's 'missing_value' attribute)
# are set to NA's.
# Argument 'signedbyte' can be TRUE for bytes to be interpreted as 
# signed, or FALSE to be unsigned.
#
ncvar_get <- function( nc, varid=NA, start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, collapse_degen=TRUE ) {

	if( class(nc) != "ncdf4" )
		stop("first argument (nc) is not of class ncdf4!")

	if( (mode(varid) != 'character') && (class(varid) != 'ncvar4') && (class(varid) != 'ncdim4') && (! is.na(varid)))
		stop(paste("Error: second argument to ncvar_get must be an object of type ncvar or ncdim",
			"(both parts of the ncdf object returned by nc_open()), the character-string name of a variable or dimension",
			"or NA to get the default variable from the file.  If the file is netcdf version 4",
			"format and uses groups, then the fully qualified var name must be given, for",
			"example, model1/run5/Temperature"))

	idobj = vobjtovarid4( nc, varid, verbose=verbose, allowdimvar=TRUE )

	have_start = (length(start)>1) || ((length(start)==1) && (!is.na(start)))
	have_count = (length(count)>1) || ((length(count)==1) && (!is.na(count)))

	#---------------------------------------------------------------------
	# Special check: if we are trying to get values from a dimvar, but the
	# dim does not have a dimvar, then just return 1:length(dim)
	#---------------------------------------------------------------------
	if( idobj$isdimvar ) {
		if( idobj$id == -1 ) {	# this happens if dim name was passed, but it has no dimvar
			#--------------------------------------------------------
			# Here we return default integers for dims with no dimvar
			#--------------------------------------------------------
			if( ! have_start )
				start <- 1
			if( ! have_count )
				count <- nc$dim[[idobj$list_index]]$len
			if( count == 1 )
				return( start )
			else
				return( start:(start+count-1) )
			}
		else
			{
			#-----------------------------------------------------------
			# Dimvars do not have list_index set, since dimvars do not
			# appear on the global var list.  However, dimvars should
			# also not have missing values, addOffsets, or scaleFactors,
			# so this is easy
			#-----------------------------------------------------------
			return( ncvar_get_inner( idobj$group_id, idobj$id, default_missval_ncdf4(), 
				start=start, count=count, verbose=verbose, signedbyte=signedbyte ))
			}
		}

	#--------------------------------------------
	# Get var's missval, addOffset, and scaleFact
	#--------------------------------------------
	if( idobj$list_index == -1 ) {
		print("internal error: list_index for var is -1!")
		print("Here is passed varid:")
		print(varid)
		}
	li = idobj$list_index
	if( nc$var[[li]]$hasAddOffset )
		addOffset = nc$var[[li]]$addOffset
	else
		addOffset = 0;
	if( nc$var[[li]]$hasScaleFact )
		scaleFact = nc$var[[li]]$scaleFact
	else
		scaleFact = 1.0;

	return( ncvar_get_inner( idobj$group_id, idobj$id, nc$var[[li]]$missval,
			addOffset, scaleFact, start=start, count=count, 
			verbose=verbose, signedbyte=signedbyte ))
}

#====================================================================================================
nc_sync <- function( nc ) {

	if( is.numeric(nc))
		ncid2use <- nc
	else if( class(nc) == 'ncdf4' )
		ncid2use <- nc$id
	else
		stop("First argument must be a simple integer ID or an object of class ncdf4, as returned by nc_open() or nc_create()")

	rv = .C("R_nc4_sync", as.integer(ncid2use), PACKAGE="pbdNCDF4")
}

#===============================================================
nc_redef <- function( nc ) {

	if( is.numeric(nc))
		ncid2use <- nc
	else if( class(nc) == 'ncdf4' )
		ncid2use <- nc$id
	else
		stop("First argument must be a simple integer ID or an object of class ncdf4, as returned by nc_open() or nc_create()")

	rv = .C("R_nc4_redef", as.integer(ncid2use), PACKAGE="pbdNCDF4")
}

#===============================================================
nc_enddef <- function( nc ) {

	if( is.numeric(nc))
		ncid2use <- nc
	else if( class(nc) == 'ncdf4' )
		ncid2use <- nc$id
	else
		stop("First argument must be a simple integer ID or an object of class ncdf4, as returned by nc_open() or nc_create()")

	rv = .C("R_nc4_enddef", as.integer(ncid2use), PACKAGE="pbdNCDF4")

	nc_sync( nc )
}

#===============================================================
nc_close <- function( nc ) {

	# Removed in version 1.5
	#if( is.numeric(nc))
	#	ncid2use <- nc

	if( class(nc) == 'ncdf4' )
		ncid2use <- nc$id
	else
		stop("First argument must be an object of class ncdf4, as returned by nc_open() or nc_create()")

	rv = .C("R_nc4_close", as.integer(ncid2use), PACKAGE="pbdNCDF4")

	#----------------------------------------------------------------------------
	# Following is taken from a posting by Simon Fear <Simon.Fear@synequanon.com>
	# to the R-help newslist on Thu, 19 Feb 2004 10:11:50 -0000
	#----------------------------------------------------------------------------
	#eval(eval(substitute(expression(nc$id <<- -1))))  # set id of CALLING object to -1
}

#===========================================================================================
# Inputs old_varname and new_varname are character strings.  If you are renaming a
# var in a group, then they both must be fully qualified varnames with the same
# number of forward slashes.
# 
# Useage:
#	ncid <- ncvar_rename( ncid, 'oldname', 'newname' )
#
ncvar_rename <- function( nc, old_varname, new_varname, verbose=FALSE ) {

	#--------------------------
	# Check inputs for validity
	#--------------------------
	if(class(nc) != "ncdf4")
		stop("Error, ncvar_rename passed something NOT of class ncdf4!")

	if( ! is.character(old_varname))
		stop("Error, ncvar_rename must be passed the old varname as a character string")

	if( ! is.character(new_varname))
		stop("Error, ncvar_rename must be passed the new varname as a character string")

	if( nslashes_ncdf4(old_varname) != nslashes_ncdf4(new_varname))
		stop("Error, if fully qualified names are passed for the old and new varnames (to rename a variable in a group), then the number of worward slashes in the old and new varnames must be the same")

	vid     <- vobjtovarid4( nc, old_varname, verbose=verbose )	# Remember, an object of class ncid4, not a simple integer
	if( vid$list_index == -1 ) 
		stop("Sorry, there was an error trying to rename the variable. Are you trying to rename a dimvar? That is not currently supported.")

	oldname <- ncvar_name( vid$group_id, vid$id )

	def_mode <- FALSE
	if( nchar(new_varname) > nchar(oldname)) {
		nc_redef( nc )	# Go back into define mode
		def_mode <- TRUE
		}

	bn_new_varname = nc4_basename( new_varname )

	rv <- list()
	rv$error <- -1
	rv <- .C("R_nc4_rename_var", 
		as.integer(vid$group_id),
		as.integer(vid$id),
		bn_new_varname,
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("error returned from C call")

	if( def_mode )
		nc_enddef( nc )

	#--------------------------------------------------------------
	# Update our netcdf class structure to reflect the changed name
	#--------------------------------------------------------------
	idx <- vid$list_index
	nc$var[[idx]]$name <- new_varname

	return(nc)
}

