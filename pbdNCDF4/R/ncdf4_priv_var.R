
#===============================================================
# Internal use only.
# Inputs:
#	ncid, varid: INTEGER values in C standard, i.e., counting starting
#		at 0, not 1.  If the var lives in a group, then ncid must
#		be the proper group ID.
# Output:
#	Character string of variable's name
#
ncvar_name <- function( ncid, varid ) {

	if( ! is.numeric(ncid))
		stop(paste("Error, must be called with numeric ncid"))

	if( ! is.numeric(varid))
		stop(paste("Error, must be called with numeric varid"))

	#------------------------------------------------------
	# This is a string NC_MAX_LEN long, to provide storage.
	# I'm not sure if this is needed or not....
	#------------------------------------------------------
	nc4_maxstring <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"
	rv <- list()
	rv$varname <- nc4_maxstring
	rv$error   <- -1
	rv <- .C("R_nc4_inq_varname",
		as.integer(ncid),
		as.integer(varid),
		varname=as.character(rv$varname),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop(paste("using ncid ",ncid," varid ",varid))

	return( rv$varname )
}

#========================================================================================
# There are numerous ways to identify a netcdf entity in the R code.  
# However the C code needs exactly one representation, an integer
# group or netcdf ID, plus an integer var ID.  This routine
# converts between the different R ways of specifiying a var
# and the actual 2 integer IDs needed by the C code. 
#
# If "allowdimvar" is true, and if the passed string does not
# match a variable's name, then we ALSO check if it matches
# a dimvar's name.  As a special check, IF the passed string matches
# a dim's name, but that dim has no dimvar, then $isdimvar is set
# to TRUE and $id is set to -1.
#
# If the input is a character string (indicating
# a variable's name) and there IS no variable of that name in 
# the file, this returns -1.
#
# Inputs:
#	nc: must be an object of class 'ncdf4'
#	varid: can be an object of class 'ncvar4' or 'ncdim4', a
#		character string, an object of class 'ncid4', or NA.
#
# Outputs:
#	an object of class 'ncid4' -- NOT A SIMPLE INTEGER!!
#	Has fields $id, $group_index, $group_id, $list_index, $isdimvar
#	(NOTE: if the input specifies a dimvar that does not actually exist in the file,
#	then $isdimvar==TRUE and $id==-1.)
# 	
vobjtovarid4 <- function( nc, varid, verbose=FALSE, allowdimvar=TRUE) {

	if( verbose )
		print(paste("vobjtovarid4: entering"))

	if( class(nc) != 'ncdf4' )
		stop('First passed argument (nc) must be an object of class ncdf4, as returned by nc_open() or nc_create()')

	if( (class(varid) != 'ncvar4') && (class(varid) != 'ncdim4') && (class(varid) != 'ncid4') && (!is.character(varid)) && (!is.na(varid)))
		stop('Second passed argument (varid) must be an object of class ncvar4 or ncid4, the character string name of a variable, or a NA (indicating to use the only var in the file)')

	#-------------------------------------------------------------
	# Easiest case is if we were given a ncid object to begin with
	#-------------------------------------------------------------
	if( class(varid) == 'ncid4' ) {
		if(verbose) print("vobjtovarid4: passed varid was an ncid, easy exit")
		return( varid ) 	# an object of class 'ncid4', not a simple integer
		}

	#------------------------------------------------------------
	# Handle case where we are given a ncvar object to begin with
	#------------------------------------------------------------
	if( class(varid) == "ncvar4" ) {

		origvarid <- varid
		if(verbose)
			print(paste("vobjtovarid4: passed a ncvar class, name=",varid$name))
		varid <- nc$var[[varid$name]]$id # Note we do NOT use varid$id in case var is from different file (but names are same)
		if( is.null(varid)) {
			print('------------------------------------------------------')
			print(paste("Error, var '", origvarid$name,"' was not found in file '", nc$filename, "'", sep=''))
			print(paste('Hint: make SURE the variable was not only defined with a call to ncvar_def(), but also included in the list of variables passed to nc_create()'))
			stop('stopping')
			}
		if( class(varid) != 'ncid4' ) {
			print('------------------------------')
			print("here is varid:")
			print(varid)
			stop(paste("Internal error #E, returned varid is not a object of class ncid4"))
			}

		#-----------------------------------------------------------------
		# Make sure this varid that will be returned has valid information
		#-----------------------------------------------------------------
		varidOK <- ((varid$id>=0) && (varid$id<=100000) && (varid$group_id >= 0))
		if( is.na(varidOK) || (!varidOK)) {
			print("vobjtovarid4: I was passed a ncvar object, BUT this object does NOT refer to any valid var in the netcdf file!")
			print(paste("This happened for netCDF filename:",nc$filename))
			print("Here are the vars in the netCDF file:")
			for( ii in 1:nc$nvars )
				print(paste(ii,": ",nc$var[[ii]]$name, sep='' ))
			print(paste("The passed varid object (which does NOT exist in that file) is:"))
			print(origvarid)
			print(paste("Hint: make SURE the variable was not only defined with a call to ncvar_def(), but also included in the list of variables passed to nc_create()"))
			stop("stopping")
			}
		
		if( verbose )
			print(paste("vobjtovarid4: returning with varid deduced from a passed ncvar object; varid=",
				varid$group_id, varid$id))
		return(varid)	# an object of class 'ncid4'
		}

	#-------------------------------------------------------------
	# Handle case where we are given an ncdim object to begin with 
	#-------------------------------------------------------------
	if( class(varid) == "ncdim4" ) {

		if( ! allowdimvar )
			stop(paste("Error, I was NOT allowed to check dimvars, but the second argument passed was an object of class ncdim4!  Name=", varid$name))

		if(verbose)
			print(paste("vobjtovarid4: passed a ncdim class, name=",varid$name))

		#-----------------------------------------------------------
		# Go through and find if there is a dim in the file that has
		# the same name as this dim.  We do not immediately use
		# this dim's dimvarid in case the dim is from a different
		# file but has the same name.
		#-----------------------------------------------------------
		name2find = varid$name
		foundit   = FALSE
		for( idim in 1:nc$ndims ) {
			if( nc$dim[[idim]]$name == name2find ) {
				#-------------------------------------------
				# Remember we return the DIMVAR, not the dim
				#-------------------------------------------
				retval  = nc$dim[[idim]]$dimvarid	# an object of type 'ncid'
				foundit = TRUE
				break
				}
			}
		if( ! foundit )
			#-----------------------------------------------------------
			# Return an ncid that indicates this is a dimvar but it does
			# not exist in the file
			#-----------------------------------------------------------
			retval = ncdf4_make_id( id=-1, group_index=-1, group_id=-1, list_index=-1, isdimvar=TRUE ) 

		if( class(retval) != 'ncid4' )
			stop(paste("Internal error #C, returned varid is not a object of class ncid4. Case with ncdim object passed; ncdim name=", varid$name))

		if( verbose )
			print(paste("vobjtovarid4: returning with varid deduced from a passed ncvar object; retval=",
				retval$group_id, retval$id))
		return(retval)	# an object of class 'ncid4'
		}

	#----------------------------------------------------------
	# If we get here, 'varid' can be a NA or a character string
	#----------------------------------------------------------

	#---------------------------------------------------------------------------
	# If varid is NA, then return the only var in the file (if there IS only one
	# var in the file).  If there is more than one var in the file, return the
	# one with the most dimensions, IF that highest-dimensionality var has more
	# dimensions than any other var in the file.  Otherwise, generate an error.
	#---------------------------------------------------------------------------
	if( (length(varid)==1) && is.na(varid)) {
		if( nc$nvars == 1 ) {
			varToUse   <- 1
			}
		else
			{
			#------------------------------------------------------------
			# Choose the most complicated var, if there is one, otherwise
			# halt with an error
			#------------------------------------------------------------
			varToUse   <- -1
			ndimsItHas <- -1
			for( ii in 1:nc$nvars ) {
				if( nc$var[[ii]]$ndims > ndimsItHas ) {
					varToUse <- ii
					ndimsItHas <- nc$var[[ii]]$ndims
					}
				}
			for( ii in 1:nc$nvars ) {
				if( (ii != varToUse) && (nc$var[[ii]]$ndims == ndimsItHas)) {
					stop(paste("File",nc$filename,"has more than one variable, so you must explicitly specify which one you want"))
					}
				}
			}
		varid <- nc$var[[varToUse]]$id	# remember, an object of class 'ncid4', not a simple int
		if( class(varid) != 'ncid4' )
			stop(paste("internal error #B, returned varid is not of class ncid4"))

		if( verbose )
			print(paste("vobjtovarid4: returning with only var in file; id=",
				nc$var[[varToUse]]$id$group_id, nc$var[[varToUse]]$id$id))
		return( varid )		# an object of class 'ncid4'
		}

	#---------------------------------------------------
	# If we get here, 'varid' must be a character string
	#---------------------------------------------------
	if( ! is.character(varid)) 
		stop("internal error: location #M: varid is not a character string!")

	origvarid <- varid

	#-------------------------------------------
	# Make sure var name follows our conventions
	#-------------------------------------------
	if( substr(varid,1,1) == '/' )
		stop(paste("Error, I was given a name that starts with a slash; fully qualified names NEVER start with a slash (this is required for backwards compatability).  Leave off the leading slash!"))

	#--------------------------------------------
	# See if any vars in this file have this name
	#--------------------------------------------
	varToUse <- -1
	if( nc$nvars > 0 ) {
		for( kk in 1:nc$nvars ) {
			if( origvarid == nc$var[[kk]]$name ) 	# check to see if fully qualified name matches
				varToUse <- kk
			}
		}

	#---------------------------------
	# Found a var with the right name,
	# return its ncid object
	#---------------------------------
	if( varToUse != -1 ) {
		if(verbose)
			print(paste("Variable named",origvarid,"found in file with varid=",
				nc$var[[varToUse]]$id$group_id, nc$var[[varToUse]]$id$id))
		varid <- nc$var[[varToUse]]$id	# remember, an object of class 'ncid4', not a simple int
		if( class(varid) != 'ncid4' ) {
			print('---- varid:')
			print(varid)
			stop(paste("internal error #A, returned varid is not of class ncid4"))
			}
		return( varid )		# an object of class 'ncid4'
		}

	#---------------------------------------------------------------
	# A var with this name was NOT found in the file.  But, it could
	# be the name of a dimvar in the file.  Check to see if we are
	# allowed to return dimvars in this case.
	#---------------------------------------------------------------
	if( ! allowdimvar ) {
		print("vobjtovarid4: error #G: I could not find the requested var in the file!")
		print(paste("requested var name:",origvarid))
		print(paste("file name:",nc$filename))
		print("Note: I was NOT allowed to check to see if this was a dimvar name")
		stop("Variable not found")
		}

	if(verbose)
		print(paste("Variable named",origvarid,"NOT found in file; looking for a dimvar with this name"))

	#-----------------------------------------------
	# Check to see if passed name matches a dim name
	#-----------------------------------------------
	for( i in 1:nc$ndims ) {
		if( origvarid == nc$dim[[i]]$name ) {
			#---------------------
			# Yes, it IS a dimvar!
			#---------------------
			varid <- nc$dim[[i]]$dimvarid 	# note: an object of class 'ncid4'.  $id will be -1 if there is no dimvar
			if( class(varid) != 'ncid4' )
				stop(paste("Internal error #D, returned varid is not a object of class ncid4"))
			if( verbose )
				print(paste("vobjtovarid4: returning with DIMvarid deduced from name; varid=",
					varid$group_id, varid$id))
			return(varid)	# an object of class 'ncid4'
			}
		}

	#------------------------------------------------------------
	# If we get here, no dimvar with the requested name was found
	#------------------------------------------------------------
	print("vobjtovarid4: error #F: I could not find the requsted var (or dimvar) in the file!")
	print(paste("var (or dimvar) name:",origvarid))
	print(paste("file name:",nc$filename))
	stop("Variable not found")
}

#===============================================================
# Internal use only.  Passed IDs must be C-style (0-based counting)
# simple integers, with the ncid being a group ID if appropriate.
# 
ncvar_inq <- function( ncid, varid ) {

	if( ! is.numeric(ncid))
		stop("Error: the first argument to ncvar_inq (ncid) is not an integer")

	if( ! is.numeric(varid))
		stop("Error: the second argument to ncvar_inq (varid) is not an integer")

	#------------------------------------------------------
	# This is a string NC_MAX_LEN long, to provide storage.
	# I'm not sure if this is needed or not....
	#------------------------------------------------------
	str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"

	rv <- list()
	rv$name    <- str.nc.max.name
	rv$varlen  <- -1
	rv$error   <- -1
	rv$type    <- -1
	rv$ndims   <- -1
	rv$natts   <- -1
	rv$precint <- -1 # INTEGER (not character) form of precision. 1=SHORT, 2=INT, 3=FLOAT, 4=DOUBLE, 5=CHAR 6=BYTE, 7=UBYTE, 8=USHORT, 9=UINT, 10=INT64, 11=UINT64, 12=STRING.  Must match C code values defined in the top of file ncdf.c
	rv$dimids  <- integer(ncvar_ndims( ncid, varid ))
	rv <- .C("R_nc4_inq_var",
		as.integer(ncid),
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
		stop("call to C function R_nc4_inq_var failed")

	var          <- list()
	var$id	     <- ncdf4_make_id( id=varid, group_index=-1, group_id=ncid, list_index=-1, isdimvar=FALSE ) 
	var$name     <- rv$name
	var$ndims    <- rv$ndims
	var$natts    <- rv$natts
	var$size     <- ncvar_size( ncid, varid )
	var$dimids   <- rv$dimids
	var$prec     <- ncvar_type_to_string( rv$precint )

	#---------------------------------------------------
	# Convert ordering of dimids from C to R conventions
	#---------------------------------------------------
	var$dimids <- var$dimids[ var$ndims:1 ]

	#--------------------------
	# Get this var's attributes
	#--------------------------
	attu <- ncatt_get_inner( ncid, varid, "units" )
	if( attu$hasatt )
		var$units <- attu$value
	else
		var$units <- ""
	attu <- ncatt_get_inner( ncid, varid, "long_name" )
	if( attu$hasatt ) 
		var$longname <- attu$value
	else
		var$longname <- var$name

	return(var)
}

#===============================================================
# Returns a vector of the size of the variable, in
# R order (XYZT).
#
# Internal use only.  Use v$varsize (where v is an object of
# class "ncvar4") if you want this info.
#
# BOTH inputs MUST be simple integers, class objects are not allowed!
# The integers are raw C-style (0-based counting), with the ncid
# actually being the group ID if necessary.
#
# A scalar variable (one with no dims) always returns a
# varsize of 1.
#
ncvar_size <- function( ncid, varid ) {

	if( mode(ncid) != 'numeric' ) 
		stop(paste("error, must be passed a numeric first arg: ncid2use, not an arg of mode", mode(ncid)))

	if( mode(varid) != 'numeric' )
		stop("Error, must be passed a numeric second arg: varid2use" )

	ndims <- ncvar_ndims( ncid, varid )
	if( ndims == 0 )
		#return(vector())	changed DWP 2012-09-20
		return(1)		# indicates a scalar var

	rv         <- list()
	rv$error   <- -1
	rv$varsize <- integer(ndims)
	rv <- .C("R_nc4_varsize",
		as.integer(ncid),
		as.integer(varid),
		varsize=as.integer(rv$varsize),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("error returned from C routine R_nc4_varsize")

	#-------------------------------------
	# Switch order from C-style to R-style
	#-------------------------------------
	rv$varsize <- rv$varsize[ ndims:1 ]

	return(rv$varsize)
}

#===============================================================
# Internal use only.   
# 
# Input: integer ncid and varid.  
#
# Output: one of the
# integer R type codes (1=short, 2=int, 3=float, 4=double,
# 5=char, 6=byte, 7=ubyte, 8=ushort, 9=uint, 10=int64, 11=uint64, 12=string).
# These are defined at the top of file ncdf.c
#
ncvar_type <- function( ncid, varid, output_string=FALSE ) {

	if( mode(ncid) != 'numeric' )
		stop("error, must be passed a numeric first arg: ncid2use")

	if( mode(varid) != 'numeric' )
		stop("Error, must be passed a numeric second arg: varid2use" )

	rv         <- list()
	rv$error   <- -1
	rv$precint <- -1

	rv <- .C("R_nc4_inq_vartype", 
		as.integer(ncid),
		as.integer(varid),	
		precint=as.integer(rv$precint),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("error returned from C call")
	return( rv$precint )
}

#======================================================================
# Internal use only.  Takes an integer precision type, converts it
# to a string.  Valus are defined at the top of file ncdf.c, and are
# not the same as in any netcdf header, since I don't want my code
# to depend on values in any netcdf library headers
#
ncvar_type_to_string = function( precint ) {

	if( precint == 1 )
		prec <- "short"
	else if( precint == 2 )
		prec <- "int"
	else if( precint == 3 )
		prec <- "float"
	else if( precint == 4 )
		prec <- "double"
	else if( precint == 5 )
		prec <- "char"
	else if( precint == 6 )
		prec <- "byte"
	else if( precint == 7 )
		prec <- "unsigned byte"
	else if( precint == 8 )
		prec <- "unsigned short"
	else if( precint == 9 )
		prec <- "unsigned int"
	else if( precint == 10 )
		prec <- "8 byte int"
	else if( precint == 11 )
		prec <- "unsinged 8 byte int"
	else if( precint == 12 )
		prec <- "string"
	else
		stop(paste("Error, unrecognized type code of variable supplied:", precint ))

	return( prec )
}

#===============================================================
# Internal use only.  Use v.ndims if you want this info.
#
# Inputs are both integers
# Output is integer number of dims in the var
#
ncvar_ndims <- function( ncid, varid ) {

	if( mode(ncid) != 'numeric' )
		stop(paste("error, ncvar_ndims must be passed a numeric first arg; mode of val passed=", mode(ncid)))

	if( mode(varid) != 'numeric' )
		stop("Error, ncvar_ndims must be passed a numeric second arg: varid2use" )

	rv <- list()
	rv$error <- -1
	rv$ndims <- -1
	rv <- .C("R_nc4_inq_varndims", 
		as.integer(ncid),
		as.integer(varid),
		ndims=as.integer(rv$ndims),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("error returned from C call")
	return( rv$ndims )
}

#===============================================================
# Internal use only
#
# Inputs: ncid = integer
#	  varname = character.   Should NOT be fully qualified
#		var name; instead the ncid should be for the proper
#		group where this var is thought to live.
#
# Returns -1 if the var is NOT found in the file, and the
# raw C-style integer varid of the var otherwise.
#
ncvar_id <- function( ncid, varname ) {

	if( mode(ncid) != 'numeric' )
		stop("error, must be passed a numeric first arg: ncid2use")

	if( mode(varname) != 'character' )
		stop("Error, must be passed a character second arg: varname" )

	if( nslashes_ncdf4(varname) > 0 )
		stop(paste("Error, var name should NOT be fully qualifed var name!!  Name=", varname))

	rv       <- list()
	rv$varid <- -1
	rv <- .C("R_nc4_inq_varid", 
		as.integer(ncid),
		as.character(varname),
		varid=as.integer(rv$varid),
		PACKAGE="pbdNCDF4")
	return(rv$varid)
}

#===========================================================================================
# This differs from the 'standard' (non-'inner') version by 
# taking ONLY the raw, C-standard correct integers for ncid
# and varid.  Note that if the variable is in a group, then
# the passed ncid must actually be the correct group ID, despite
# its name.  
#
# 'missval' is the value of the variable's missing value attribute
# (or NA if the variable does not have such an attribute).  Any
# read-in values that are equal to missval are set to NA.
#
ncvar_get_inner <- function( ncid, varid, missval, addOffset=0., scaleFact=1.0, start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, collapse_degen=TRUE ) {

	if( ! is.numeric(ncid))
		stop("Error, first arg passed to ncvar_get_inner (ncid) must be a simple C-style integer that is passed directly to the C api")
	if( ! is.numeric(varid))

		stop("Error, second arg passed to ncvar_get_inner (varid) must be a simple C-style integer that is passed directly to the C api")

	if( verbose ) 
		print(paste("ncvar_get_inner: entering with (C-STYLE INTEGER ONLY) ncid=", ncid, 
			"varid=", varid ))

	tmp_typename = c('short', 'int', 'float', 'double', 'char', 'byte' )

	have_start = (length(start)>1) || ((length(start)==1) && (!is.na(start)))
	have_count = (length(count)>1) || ((length(count)==1) && (!is.na(count)))

	sm <- storage.mode(start)
	if( (sm != "double") && (sm != "integer") && (sm != "logical"))
		stop(paste("passed a start argument of storage mode",sm,"; can only handle double or integer"))
	sm <- storage.mode(count)
	if( (sm != "double") && (sm != "integer") && (sm != "logical"))
		stop(paste("passed a 'count' argument with storage mode '",sm,"'; can only handle double or integer", sep=''))

	if( signedbyte )
		byte_style = 1	# 1=signed
	else
		byte_style = 2	# 2=unsigned

	varsize <- ncvar_size ( ncid, varid )
	ndims   <- ncvar_ndims( ncid, varid )
	if( verbose ) {
		print(paste("ndims:",ndims))
		print("ncvar_get: varsize:")
		print(varsize)
		}

	#------------------------------
	# Fix up start and count to use
	#------------------------------
	if( ndims == 0 ) {
		start <- 1
		count <- 1
		}
	else
		{
		if( ! have_start )
			start <- rep(1,ndims)	# Note: use R convention for now
		if( ! have_count )
			count <- varsize - start + 1	
		else
			{
			#------------------
			# Take care of -1's
			#------------------
			count <- ifelse( (count == -1), varsize-start+1, count)
			}
		}
	if( verbose ) {
		print("ncvar_get: start:")
		print(start)
		print("ncvar_get: count:")
		print(count)
		}

	if( ndims > 0 ) {
		if( length(start) != ndims ) 
			stop(paste("Error: variable has",ndims,"dims, but start has",length(start),"entries.  They must match!"))
		if( length(count) != ndims ) 
			stop(paste("Error: variable has",ndims,"dims, but count has",length(count),"entries.  They must match!"))
		}

	#----------------------------------------
	# Need to know how much space to allocate
	#----------------------------------------
	totvarsize <- prod(count)
	if( verbose )
		print(paste("ncvar_get: totvarsize:",totvarsize))
	
	#--------------------------------------------------
	# Switch from R to C convention for start and count
	#--------------------------------------------------
	c.start <- start[ ndims:1 ] - 1
	c.count <- count[ ndims:1 ]

	rv <- list()
	rv$error <- -1

	#---------------------------------
	# Get the correct type of variable
	#---------------------------------
	precint <- ncvar_type( ncid, varid ) # 1=short, 2=int, 3=float, 4=double, 5=char, 6=byte, 7=ubyte, 8=ushort, 9=uint, 10=int64, 11=uint64, 12=string
	if( verbose )
		print(paste("Getting var of type",tmp_typename[precint]))
	if( (precint == 1) || (precint == 2) || (precint == 6) || (precint == 7) || (precint == 8) || (precint == 9)) {
		#--------------------------------------
		# Short, Int, Byte, UByte, UShort, Uint
		#--------------------------------------
		rv$data  <- integer(totvarsize)
		rv <- .C("R_nc4_get_vara_int", 
			as.integer(ncid),
			as.integer(varid),	
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			as.integer(byte_style), # 1=signed, 2=unsigned
			data=as.integer(rv$data),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4",
			DUP=TRUE)
		if( rv$error != 0 ) 
			stop("C function R_nc4_get_var_int returned error")
		}
	else if( (precint == 3) || (precint == 4) || (precint == 10) || (precint == 11)) {
		#-----------------------------
		# Float, double, int64, uint64
		#-----------------------------
		rv$data  <- double(totvarsize)
		rv <- .C("R_nc4_get_vara_double", 
			as.integer(ncid),
			as.integer(varid),
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.double(rv$data),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4",
			DUP=TRUE)
		if( rv$error != 0 ) 
			stop("C function R_nc4_get_vara_double returned error")
		}
	else if( precint == 5 ) {
		#-----
		# Char 
		#-----
		strndims <- ndims - 1
		strlen   <- count[1] + 1
		strdim   <- 1
		if( strndims >= 1 ) {
			strdim <- count[2:ndims]
			nstr   <- prod(strdim)
			}
		else
			nstr <- 1
		if(verbose)
			print(paste("ndims:",ndims,"strndims:",strndims,"strlen:",strlen,"nstr:",nstr))

		#----------------------------------------------
		# Make a character string of the specified size
		#----------------------------------------------
		stor     <- blankstring_ncdf4( totvarsize )
		stordata <- blankstring_ncdf4(strlen)
		if( verbose )
			print(paste("length of stor string:",nchar(stor)))
		rv$tempstore <- stor
		rv$data      <- array(stordata, dim=strdim)

		rv <- .C("R_nc4_get_vara_text", 
			as.integer(ncid),
			as.integer(varid),
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			tempstore=as.character(rv$tempstore),
			data=as.character(rv$data),
			error=as.integer(rv$error),
			PACKAGE="pbdNCDF4")
		if( rv$error != 0 ) 
			stop("C function R_nc4_get_var_text returned error")

		dim(rv$data) <- strdim
		}

	else if( precint == 12 ) {
		#-----------------------------
		# netcdf version 4 String type
		#-----------------------------
		rv <- .Call( "R_nc4_get_vara_string",
			as.integer(ncid),
			as.integer(varid),
			as.integer(c.start),    # Already switched to C convention...
			as.integer(c.count),    # Already switched to C convention...
			PACKAGE="pbdNCDF4" )
		}
	else
		{
		stop(paste("Trying to get variable of an unhandled type code: ",precint, "(", ncvar_type_to_string(precint), ")"))
		}
	if( verbose )
		print(paste("ncvar_get: C call returned",rv$error))

	#--------------------------------------------------------
	# Set our dims...but collapse degenerate dimensions first
	#--------------------------------------------------------
	if( ndims > 0 ) {
		if( collapse_degen ) {
			count.nodegen <- vector()
			foundone <- 0
			for( i in 1:ndims )
				if( count[i] > 1 ) {
					count.nodegen <- append(count.nodegen, count[i])
					foundone <- 1
					}
			if( foundone == 0 ) 
				dim(rv$data) <- (1)
			else
				{
				if( verbose )
					print(paste("count.nodegen:",count.nodegen,"   Length of data:",length(rv$data)))
				if( precint != 5 )
					dim(rv$data) <- count.nodegen
				}
			}
		if( verbose ) {
			print("ncvar_get: final dims of returned array:")
			print(dim(rv$data))
			}
		}

	#----------------------------------------------------------
	# Change missing values to "NA"s.  Note that 'varid2Rindex'
	# is NOT filled out for dimvars, so skip this if a dimvar
	# 1=short, 2=int, 3=float, 4=double, 5=char, 6=byte
	#----------------------------------------------------------
	if( precint != 5 ) {
		if( verbose ) print("ncvar_get: setting missing values to NA")
		if( (precint==1) || (precint==2) || (precint==6) || (precint==7) || (precint==8) || (precint==9)) {
			#--------------------------------------
			# Short, Int, Byte, UByte, UShort, UInt
			#--------------------------------------
			if( verbose ) print(paste("ncvar_get_inner: setting ", tmp_typename[precint],"-type missing value of ", missval, " to NA", sep=''))
			if( ! is.na(missval) ) 
				rv$data[rv$data==missval] <- NA
			}
		else if( (precint==3) || (precint==4) || (precint==10) || (precint==11)) {
			#-----------------------------------------------
			# Float, Double, 8-byte int, unsigned 8-byte int
			#-----------------------------------------------
			if( ! is.na(missval) ) {
				tol <- abs(missval*1.e-5)
				if( verbose ) print(paste("ncvar_get_inner: setting ", tmp_typename[precint],"-type missing value of ", missval, 
					" (w/tolerance ", tol,") to NA", sep=''))
				rv$data[abs(rv$data-missval)<tol] <- NA
				}
			}
		}

	#--------------------------------------
	# Implement add_offset and scale_factor
	#--------------------------------------
	if( (scaleFact != 1.0) || (addOffset != 0.0) ) {
		if( verbose ) 
			print(paste("ncvar_get: implementing add_offset=", addOffset, " and scaleFact=", scaleFact ))
		rv$data <- rv$data * scaleFact + addOffset
		}

	return(rv$data)
}

#=======================================================================================================
ncvar_def_deflate = function( root_id, varid, shuffle, deflate, deflate_level ) {

	if( !is.numeric(root_id))
		stop("must be passed a numeric root_id")
	if( !is.numeric(varid))
		stop("must be passed a numeric varid")
	if( !is.numeric(shuffle))
		stop("must be passed a numeric shuffle")
	if( !is.numeric(deflate))
		stop("must be passed a numeric deflate")
	if( !is.numeric(deflate_level))
		stop("must be passed a numeric deflate_level")

	shuffle = as.integer(shuffle)
	if( (shuffle != 0) && (shuffle != 1))
		stop("shuffle must be either 0 or 1")

	deflate = as.integer(deflate)
	if( (deflate != 0) && (deflate != 1))
		stop("deflate must be either 0 or 1")

	deflate_level = as.integer(deflate_level)
	if( (deflate < 0) || (deflate > 9))
		stop("deflate_level must be between 0 and 9, inclusive")

	rv         <- list()
	rv$root_id <- root_id
	rv$varid   <- varid
	rv$shuffle <- shuffle
	rv$deflate <- deflate
	rv$deflate_level <- deflate_level
	rv$error   <- -1
	rv <- .C("R_nc4_def_var_deflate", 
		as.integer(rv$root_id),
		as.integer(rv$varid),
		as.integer(rv$shuffle),
		as.integer(rv$deflate),
		as.integer(rv$deflate_level),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("C function R_nc4_def_var_deflate returned error")
}

#=======================================================================================================
# Returns a list with:
#	$shuffle: ==1 if shuffle is turned on for this var, ==0 otherwise
#	$deflate: ==1 if deflate (compression) filter is turned on for this var, ==0 otherwise
#	$deflate_level: if $deflate==1, then this integer (0-9) indicates the deflate (compression) level
#
# Note that this routine must ONLY be called for variables in a NETCDF-4 format file!  
#
ncvar_inq_deflate = function( root_id, varid ) {

	if( !is.numeric(root_id))
		stop("must be passed a numeric root_id")
	if( !is.numeric(varid))
		stop("must be passed a numeric varid")

	rv         <- list()
	rv$root_id <- root_id
	rv$varid   <- varid
	rv$shuffle <- 0
	rv$deflate <- 0
	rv$deflate_level <- 0
	rv$error   <- -1
	rv <- .C("R_nc4_inq_var_deflate", 
		as.integer(rv$root_id),
		as.integer(rv$varid),
		shuffle=as.integer(rv$shuffle),
		deflate=as.integer(rv$deflate),
		deflate_level=as.integer(rv$deflate_level),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("C function R_nc4_inq_var_deflate returned error")

	retval = list( shuffle=rv$shuffle, deflate=rv$deflate, deflate_level=rv$deflate_level )

	return( retval )
}

#=======================================================================================================
# NOTE: on entry, length(chunksizes) MUST EQUAL ndims in the var.  This routine assumes
# this is true and will crash otherwise.  So check this before calling this routine!
#
# 'chunksizes' array is in R order upon entry.  In this routine we convert it to C order.
#
ncvar_def_chunking = function( root_id, varid, storage, chunksizes ) {

	if( !is.numeric(root_id))
		stop("must be passed a numeric root_id")
	if( !is.numeric(varid))
		stop("must be passed a numeric varid")
	if( !is.numeric(storage))
		stop("must be passed a numeric storage")
	if( !is.numeric(chunksizes))
		stop("must be passed a numeric chunksizes")

	chunksizes = as.integer(chunksizes)
	chunksizes = chunksizes[length(chunksizes):1]	# Switch from R order to C order

	storage = as.integer(storage)
	if( (storage != 1) && (storage != 2))
		stop("storage must be either 1 (for NC_CONTIGUOUS) or 2 (for NC_CHUNKED)")

	rv         <- list()
	rv$root_id <- root_id
	rv$varid   <- varid
	rv$ndims   <- length(chunksizes)
	rv$storage <- storage
	rv$chunksizes <- chunksizes
	rv$error   <- -1
	rv <- .C("R_nc4_def_var_chunking", 
		as.integer(rv$root_id),
		as.integer(rv$varid),
		as.integer(rv$ndims),
		as.integer(rv$storage),
		as.integer(rv$chunksizes),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("C function R_nc4_def_var_chunking returned error")
}

#=======================================================================================================
# Returns a list with $storage and $chunksizes.  
# $storage is 1 for contiguous storage and 2 for chunked storage.
# $chunksizes is an array of length ndims, in R order.
#
ncvar_inq_chunking = function( root_id, varid, ndims ) {

	if( !is.numeric(root_id))
		stop("must be passed a numeric root_id")
	if( !is.numeric(varid))
		stop("must be passed a numeric varid")
	if( !is.numeric(ndims))
		stop("must be passed a numeric storage")

	chunksizes = array( as.integer(0), ndims )

	rv         <- list()
	rv$root_id <- root_id
	rv$varid   <- varid
	rv$ndims   <- ndims
	rv$storage <- 0
	rv$chunksizes <- array(as.integer(0), ndims)
	rv$error   <- -1
	rv <- .C("R_nc4_inq_var_chunking", 
		as.integer(rv$root_id),
		as.integer(rv$varid),
		as.integer(rv$ndims),
		storage=as.integer(rv$storage),
		chunksizes=as.integer(rv$chunksizes),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("C function R_nc4_inq_var_chunking returned error")

	#-------------------------------------------------------------
	# NOTE we switch chunksizes from C ordering to R ordering here
	# (i.e., we reverse it)
	#-------------------------------------------------------------
	retval = list( storage=rv$storage, chunksizes=rv$chunksizes[ndims:1] )

	return( retval )
}

