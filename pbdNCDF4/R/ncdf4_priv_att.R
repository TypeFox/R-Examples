
#===============================================================
# FOR INTERNAL USE ONLY.
# ASSUMES nc and varid are simple C-style integers.
# Set varid to -1 to get the global # of atts.
#
ncatt_get_n <- function( nc, varid ) {

	verbose = FALSE

	if( verbose ) print(paste("ncatt_get_n: entering with integer (ONLY) nc=", nc, "and integer (ONLY) varid=",varid))

	if( ! is.numeric(nc))
		stop(paste("error, first arg must be of class ncdf4!"))

	if( ! is.numeric(varid))
		stop(paste("Error, second arg must be an integer!"))

	if( varid == -1 ) {
		if( verbose ) print ("ncatt_get_n: varid == -1, so getting number of global attributes")
		rv            <- list()
		rv$ndims      <- -1
		rv$nvars      <- -1
		rv$natts      <- -1
		rv$error      <- -1
		rv <- .C("R_nc4_inq",
			as.integer(nc),
			ndims=as.integer(rv$ndims),
			nvars=as.integer(rv$nvars),
			natts=as.integer(rv$natts),
			## REMOVE unlimdimid=as.integer(rv$unlimdimid),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		if( rv$error != 0 ) 
			stop(paste("R_nc4_inq returned error on file",nc$filename,"!"))
		}
	else
		{
		if( verbose ) print ("ncatt_get_n: varid != -1, so getting number of attributes for a specific var")
		str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"
		rv <- list()
		rv$name    <- str.nc.max.name
		rv$varlen  <- -1
		rv$error   <- -1
		rv$type    <- -1
		rv$ndims   <- -1
		rv$natts   <- -1
		rv$precint <- -1 
		rv$dimids  <- integer(ncvar_ndims( nc, varid ))
		rv <- .C("R_nc4_inq_var",
			as.integer(nc),
			as.integer(varid),
			name=as.character(rv$name),
			type=as.integer(rv$type),
			ndims=as.integer(rv$ndims),
			dimids=as.integer(rv$dimids),
			natts=as.integer(rv$natts),
			precint=as.integer(rv$precint),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		if( rv$error != 0 ) 
			stop(paste("R_nc4_inq_var returned error on file",nc$filename,"!"))
		}

	return( rv$natts )
}

#======================================================================================================
# Args ncid and varid passed to this routine must be simple integers that can be directly
# passed to the C API.  I.e., 0-based counting, and a global attribute is flagged by -1.
# NOTE: putting a global attribute is flagged by varid == -1 (NOT zero, since args passed to this
# function must be the C API values, not the R values!)
#
ncatt_put_inner = function( ncid, varid, attname, attval, prec=NA, verbose=FALSE, definemode=FALSE ) {

	if( verbose ) print(paste('ncatt_put_inner: entering with ncid=', ncid, 'varid=', varid, 'attname=', attname, 'attval=', attval, 'prec=', prec, 'storage.mode(attval)=', storage.mode(attval) ))

	if( ! is.numeric(ncid))
		stop("Can only be called with a simple C-style (0-based counting) integer ncid")

	if( ! is.numeric(varid))
		stop("Can only be called with a simple C-style (0-based counting) integer varid")

	if( is.null(attval)) {
		print(paste("Warning: ncatt_put passed a NULL attribute; name=", attname ))
		return()
		}

	#-------------------------------------------------------------
	# Note there are TWO types here.  One is the storage mode
	# of the passed attval.  The netCDF routine to call is based
	# on this stoarge mode. The second type is the type of 
	# attribute to create.  This is passed as the 'prec' parameter 
	# to this netcDF routine.  If the passed prec is NA, then
	# the storage mode of the attval is used as the type of the
	# attribute to create.
	#-------------------------------------------------------------

	#---------------------------------------------------------
	# Get the netCDF function to call ... this always depends
	# exclusively on the storage mode of the attval
	#---------------------------------------------------------
	if( storage.mode(attval) == "integer" ) 
		funcname <- "R_nc4_put_att_int"
	else if( storage.mode(attval) == "double" )
		funcname <- "R_nc4_put_att_double"
	else if( storage.mode(attval) == "character")
		funcname <- "R_nc4_put_att_text"
	else if( storage.mode(attval) == "logical")
		funcname <- "R_nc4_put_att_logical"
	else
		stop(paste("ncatt_put: error, passed an attribute with a storage mode not handled.  Att name:",attname,"Att value:",attval,"Storage mode passed:",storage.mode(attval),".  Handled types: integer double character"))

	if( verbose ) print(paste("ncatt_put_inner: using function",funcname))

	#-------------------------------------------------------------
	# Get the type of attribute to create.  This follows the var's
	# type, in general, but can be manually overridden.
	#-------------------------------------------------------------
	atttypeShort <- 1  # These MUST match the values in the C code
	atttypeInt   <- 2 
	atttypeFloat <- 3
	atttypeDbl   <- 4
	atttypeText  <- 5
	atttypeByte  <- 6
	typetocreate <- -1
	if( (length(prec)==1) && is.na(prec) ) {
		if( verbose ) print(paste("ncatt_put_inner: no user-specified att type was given, figuring it out..."))

		#---------------------------------------------------------------------------
		# The logic of this code is as follows. In general, given no additional
		# information, it would be nice for the attribute to be the same type as
		# the original variable IF POSSIBLE.  Now, if we are given a global 
		# attribute, there is no "original variable", so we just take the precision
		# (type) of the attribute as it is. If there IS a variable, we would like
		# to make the att the same type if they are "compatible" types. For example,
		# if the var is an int and the passed att is 52.0, it would be nice to 
		# store it as an int rather than as a float. This code attempts to 
		# make that decision.  In general, this only applies if we are trying to
		# cast near-integer floats or doubles to int in the case of an int var.
		#---------------------------------------------------------------------------
		if( varid == -1 ) 	# A global attribute
			prec <- storage.mode(attval)

		else if( storage.mode(attval) == "character" ) {
			prec <- "character"	# This always works, but MAY be inconvenient. However the user passed a char, so honor that request!
			}

		else
			{
			prec = storage.mode(attval)	# our default choice

			#----------------------------------------------------------------
			# Get the prec (type) of the VARIABLE this att is associated with
			#----------------------------------------------------------------
			var_precint = ncvar_type( ncid, varid )
			var_prec    = ncvar_type_to_string( var_precint )

			if( var_prec == "int" ) {
				att_is_int = (is.numeric(attval) && (floor(attval) == attval))
				if( att_is_int ) 
					prec = 'int'
				}
			}
		if( verbose ) print(paste("ncatt_put_inner: using deduced attribute prec of", prec))
		}
	else
		if( verbose ) print(paste("ncatt_put_inner: using specified attribute prec of", prec))

	if( verbose ) print(paste("ncatt_put_inner: prec to create:",prec))
	if( (prec == "single") || (prec == "float"))
		typetocreate <- atttypeFloat
	else if( prec == "short" )
		typetocreate <- atttypeShort
	else if( prec == "byte" )
		typetocreate <- atttypeByte
	else if( prec == "double" )
		typetocreate <- atttypeDbl
	else if( (prec == "integer" ) || (prec == "int"))
		typetocreate <- atttypeInt
	else if( (prec == "text") || (prec == "character") || (prec == "char"))
		typetocreate <- atttypeText
	else
		stop(paste("Error in ncatt_put: unknown prec type specified:",prec,". Known values: short integer float double character"))

	if( ! definemode )
		nc_redef(ncid)

	rv     <- list()
	rv$error <- -1
	### WCC: R CMD check warning.
	# rv <- .C(funcname,
	# 	as.integer(ncid),
	# 	as.integer(varid),
	# 	as.character(attname), 
	# 	as.integer(typetocreate),
	# 	as.integer(length(attval)),
	# 	attval,
	# 	error=as.integer(rv$error),
	# 	PACKAGE="pbdNCDF4",
	# 	NAOK=TRUE )
	### WCC: avoid R CMD check warning.
	if( storage.mode(attval) == "integer" ){
		rv <- .C("R_nc4_put_att_int",
	 		as.integer(ncid),
	 		as.integer(varid),
	 		as.character(attname), 
	 		as.integer(typetocreate),
	 		as.integer(length(attval)),
	 		attval,
	 		error=as.integer(rv$error),
	 		PACKAGE="pbdNCDF4",
	 		NAOK=TRUE )
	} else if( storage.mode(attval) == "double" ){
		rv <- .C("R_nc4_put_att_double",
	 		as.integer(ncid),
	 		as.integer(varid),
	 		as.character(attname), 
	 		as.integer(typetocreate),
	 		as.integer(length(attval)),
	 		attval,
	 		error=as.integer(rv$error),
	 		PACKAGE="pbdNCDF4",
	 		NAOK=TRUE )
	} else if( storage.mode(attval) == "character" ){
		rv <- .C("R_nc4_put_att_text",
	 		as.integer(ncid),
	 		as.integer(varid),
	 		as.character(attname), 
	 		as.integer(typetocreate),
	 		as.integer(length(attval)),
	 		attval,
	 		error=as.integer(rv$error),
	 		PACKAGE="pbdNCDF4",
	 		NAOK=TRUE )
	} else if( storage.mode(attval) == "logical" ){
		rv <- .C("R_nc4_put_att_logical",
	 		as.integer(ncid),
	 		as.integer(varid),
	 		as.character(attname), 
	 		as.integer(typetocreate),
	 		as.integer(length(attval)),
	 		attval,
	 		error=as.integer(rv$error),
	 		PACKAGE="pbdNCDF4",
	 		NAOK=TRUE )
	} else{
		stop(paste("ncatt_put: error, passed an attribute with a storage mode not handled.  Att name:",attname,"Att value:",attval,"Storage mode passed:",storage.mode(attval),".  Handled types: integer double character"))
	}

	if( rv$error != 0 ) {
		print(paste("Error in ncatt_put, while writing attribute",
			attname,"with value",attval))
		stop(paste("Error return from C call",funcname,"for attribute",attname))
		}

	if( ! definemode )
		nc_enddef(ncid)
}

#===================================================================================================
# The difference between the "inner" version and the "regular" version is that the inner
# version is passed only simple C-style integer ID's to operate on.  
# Inputs:
#	ncid, varid: C-style (0-based counting) integers.  ncid must actually be a group
#		ID if appropriate.
#
# SPECIAL NOTE:: Ordinarily, the C interface uses a varid == -1 to indicate global attributes.
# while the R code visible to the user indicates global attributes using a varid == 0.  Since this
# routine takes as its input the actual C-style numbers, the passed varid must equal -1 to access 
# global attributes.
#
ncatt_get_inner <- function( ncid, varid, attname=NA, verbose=FALSE ) {

	if( ! is.numeric(ncid))
		stop(paste("ncatt_get_inner must be passed a simple C-style (0-based counting) integer as the first argument (ncid)"))

	if( ! is.numeric(varid))
		stop(paste("ncatt_get_inner must be passed a simple C-style (0-based counting) integer as the first argument (ncid)"))
	
	retval <- list()
	
	#---------------------------------------------------------------------------
	# If no attname is specified, return a list with attribute name/value pairs
	#---------------------------------------------------------------------------
	if( is.na(attname)) {
		if( verbose ) print(paste("ncatt_get_inner: no attname specified, returning a list with name/value pairs *******"))
		na <- ncatt_get_n( ncid, varid )
		if( verbose ) print(paste("ncatt_get_inner: number of atts for this var [or file, if global]:", na))
		if( na == 0 ) {
			if( verbose ) print(paste("ncatt_get_inner: no attributes for this var/file, returning empty list"))
			return( retval )
			}
		if( verbose ) print(paste("ncatt_get_inner: Looping over",na,"attributes ******"))
		for( iatt in 0:(na-1) ) {	# NOTE C-style 0-based counting here
			#---------------------
			# Get attribute's name
			#---------------------
			str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"
			rv0 <- list()
			rv0$error  <- -1
			rv0$attname <- str.nc.max.name
			rv0 <- .C("R_nc4_inq_attname",
				as.integer(ncid),
				as.integer(varid),
				as.integer(iatt),	
				attname=as.character(rv0$attname),
				error=as.integer(rv0$error),
				PACKAGE="pbdNCDF4")
			if( rv0$error != 0 ) {
				stop("error on call to R_nc4_inq_attname")
				}
			if( verbose ) print(paste("ncatt_get_inner: Attribute", iatt, "(in 0-based counting) has name \"", rv0$attname, "\"" ))

			#-------------------------------------------------------------------------------
			# Get attribute's value by recursively calling myself with the name now supplied
			#-------------------------------------------------------------------------------
			if( verbose ) print(paste("ncatt_get_inner: recursively calling myself to value for attribute \"", rv0$attname, "\"" ))
			tt <- ncatt_get_inner( ncid, varid, rv0$attname )
			if( ! tt$hasatt ) 
				stop(paste("internal error: could not get attribute value for att named ", rv0$attname ))
			retval[[rv0$attname]] <- tt$value
			}
		if( verbose ) print(paste("ncatt_get_inner: done, returning a list with name/value pairs *******"))
		return( retval )
		}

	#----------------------------------------------------
	# Find out if the attribute exists for this variable, 
	# and, if so, what type and length it is.
	#----------------------------------------------------
	rv0 <- list()
	rv0$error  <- -1
	rv0$attlen <- -1
	rv0$type   <- -1
	rv0 <- .C("R_nc4_inq_att",
		as.integer(ncid),
		as.integer(varid),
		as.character(attname),
		type=as.integer(rv0$type), # 1=short 2=int 3=float 4=double 5=text 6=byte 7=ubyte 8=ushort 9=uint 10=int64 11=uint64 12=string
		attlen=as.integer(rv0$attlen),
		error=as.integer(rv0$error),
		PACKAGE="pbdNCDF4")
	if( rv0$error != 0 ) {
		#---------------------------------------------------------
		# This variable did NOT have an attribute named 'attname',
		# or it is of a type not handled.
		#---------------------------------------------------------
		retval$hasatt <- FALSE
		retval$value  <- 0
		return(retval)
		}

	retval$hasatt <- TRUE

	rv <- list()
	rv$error     <- -1

	if( (rv0$type == 1) || (rv0$type == 2) || (rv0$type == 6) || (rv0$type == 7) || (rv0$type == 8) || (rv0$type == 9)) {
		#--------------------------------------
		# Short, Int, Byte, UByte, UShort, UInt 
		#--------------------------------------
		rv$attribute <- rep(as.integer(0),rv0$attlen)
		rv <- .C("R_nc4_get_att_int",
			as.integer(ncid),
			as.integer(varid),
			as.character(attname),
			attribute=as.integer(rv$attribute),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		}
	else if( (rv0$type == 3) || (rv0$type == 4) || (rv0$type == 10) || (rv0$type == 11)) {
		if( (rv0$type == 10) || (rv0$type == 11)) {
			print(paste(">>>> WARNING <<<  attribute", attname, "is an 8-byte value, but R"))
			print(paste("does not support this data type. I am returning a double precision"))
			print(paste("floating point, but you must be aware that this could lose precision!"))
			}
		#------------------------------------------------
		# Single, Double, 8-byte int, unsigned 8-byte int
		#------------------------------------------------
		rv$attribute <- rep(0.0,rv0$attlen)
		rv <- .C("R_nc4_get_att_double",
			as.integer(ncid),
			as.integer(varid),
			as.character(attname),
			attribute=as.double(rv$attribute),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		}
	else if( rv0$type == 5 ) {
		#-----------------------------------------------------------
		# Character string ... note we allocate storage for it first
		#-----------------------------------------------------------
		rv$attribute <- blankstring_ncdf4( rv0$attlen )
		rv <- .C("R_nc4_get_att_text",
			as.integer(ncid),
			as.integer(varid),
			as.character(attname),
			attribute=as.character(rv$attribute),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		}
	else
		stop("error, unhandled attribute type!")

	if( rv$error != 0 ) {
		#---------------------------------------------
		# ? Got some strange error -- return as if the
		# attribute did not exist
		#---------------------------------------------
		retval$hasatt <- FALSE
		retval$value  <- 0
		return(retval)
		}

	retval$value <- rv$attribute
	if( verbose ) print(paste("ncatt_get_inner: done for a single attribute, value= >",retval$value,"<"))
	return(retval)
}


