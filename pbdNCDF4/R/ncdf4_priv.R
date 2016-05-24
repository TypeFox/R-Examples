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
# 9-April-2001
#
#-----------------------------------------------------------------
# Netcdf-version 4 notes.
# Not much needs to be changed in the interface to make it work
# for netcdf version 4.  Basically, if you have a slash in the
# variable or dim name, such as "models/Temperature", then "models"
# is a group that will be auto-created as necessary.  Same goes
# for dims.  
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
#	varid2Rindex : for internal use only; maps a numeric varid stright
#		from the netcdf file to which element in the list of vars this var is.
#	writable: TRUE or FALSE
#
# class: ncdim4 (returned by dim.def.ncdf, which creates a NEW 
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
# class: ncvar4 (returned by var.def.ncdf, which creates a NEW
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
# class: vals.ncdf (returned by ncvar_get)
#	xvals, yvals, zvals, tvals: the dimensional values, as appropriate
#	vals	: the data values
#

#=================================================================
# Needed because the simpler construction 
#
#	for( i in 1:n )
#
# fails if n==0.  Instead, use
#
#	for( i in nc4_loop(1,n) ) 
#
# if there is any possibility that n can be zero.
#
nc4_loop <- function( a, b ) {

	if( b < a )
		return(NULL);
	return( a:b );
}

#=================================================================
# Given a name of the form "dir/dir/file", or
# "/dir/dir/file", or "file", or "/file", this
# returns the dirname (stuff before the file) or
# basename (filename).  This version does not 
# strip suffixes from filenames, as I don't need
# that functionality.
#
# Input string		dir=FALSE	dir=TRUE
# ------------		---------	---------
# /dir1/fname		fname		/dir1
# /dir1/dir2/fname	fname		/dir1/dir2
# dir1/fname		fname		dir1
# dir1/dir2/fname	fname		dir1/dir2
# /fname		fname		"" (empty string)
# fname			fname		"" (empty string)
#
nc4_basename <- function( nam, dir=FALSE ) {

	if( ! dir ) {
		ss <- strsplit( nam, '/' )
		ns <- length( ss[[1]] )
		return( ss[[1]][ns] )
		}

	#------------------------------
	# If get here, doing 'dir=TRUE'
	#------------------------------
	nmax     <- nchar( nam )
	nslash   <- 0
	slashloc <- array(0,nmax)
	for( i in 1:nchar(nam)) {
		if( substr( nam, i, i ) == '/' ) {
			nslash <- nslash + 1
			slashloc[nslash] <- i
			}
		}

	if( nslash == 0 )
		return( "" )

	return( substr( nam, 1, slashloc[nslash]-1 ))
}

#=================================================================
# Utility to return a streng of length 'n'; this is used for storage
#
blankstring_ncdf4 <- function( n ) {

	rv <- .Call("R_nc4_blankstring", as.integer(n), PACKAGE="pbdNCDF4")
	return( rv )
}

#=================================================================
# This is where the default missing value for vars is set.
#
default_missval_ncdf4 <- function()  
{
	return(1.e30)
}

#=================================================================
# Returns the number of forward slash characters ("/") in a name
#
nslashes_ncdf4 <- function( nam ) {

	nfs <- 0	# number of forward slashes

	for( ic in nc4_loop(1,nchar(nam)) )
		if( substr(nam,ic,ic) == '/' )
			nfs <- nfs + 1

	return( nfs )
}

#===============================================================
# In netcdf v4, there can be multiple vars with the 
# same name and same ID, if they are in different groups.  So 
# we have to have a more complex object for the var "ID" than just 
# a simple integer.  Note that this does not seem to be true for
# dims, which, as far as I can tell, are numbered consecutively 
# throughout a netcdf file that has groups.
#
# Elements:
#	id: the integer varid or dimid to use when calling the 
#		C routines to access the dim/var WITH THE CORRECT GROUP ID
#	group_id: the group ID to use when accessing the C routines to
#		manipulate this dim or var
#	group_index: the index in to the nc$group list for the group where
#		this dim or var lives
#	list_index: this is a number N, between 1 and nc$nvars (if for a
#		variable) or between 1 and nc$ndims (if for a dim)
#		that shows the entry of this dim/var on the overall
#		ncdf $var or $dim list.  For example, if list_index==3
#		for a variable, then nc$var[[3]] is the variable we
#		are talking about.
#	isdimvar: TRUE if a dimvar, FALSE otherwise
#
# Note that none of these is valid until the dim or var is 
# actually created on disk; this is indicted by the default vals
# of -1.
#
ncdf4_make_id <- function( id=-1, group_index=-1, group_id=-1, list_index=-1, isdimvar=FALSE ) {

	retval	  	  <- list( id=id, group_index=group_index, group_id=group_id,
				list_index=list_index, isdimvar=isdimvar )
	class(retval)  	  <- 'ncid4'
	return( retval )
}

#======================================================================================
# Given a netcdf ID, returns the netcdf format.  Values here MUST MATCH THE VALUES
# COMPILED INTO THE C CODE, file ncdf.c, routine R_nc4_inq_format, about line 1123
# or so.
#
ncdf4_format <- function( root_id ) {

	if( ! is.integer(root_id))
		root_id <- as.integer( root_id + 0.1 )

	ierr <- as.integer(-1)
	rv <- .Call( "R_nc4_inq_format", as.integer(root_id), ierr, PACKAGE="pbdNCDF4" )
	if( ierr != 0 ) {
		stop(paste("Error in nc_grpname, encountered when root_id=",root_id))
		}

	if( rv == 1 )
		return( 'NC_FORMAT_CLASSIC' )
	else if( rv == 2 )
		return( 'NC_FORMAT_64BIT' )
	else if( rv == 3 ) 
		return( 'NC_FORMAT_NETCDF4' )
	else if( rv == 4 ) 
		return( 'NC_FORMAT_NETCDF4_CLASSIC' )
	else
		stop(paste("C call R_nc4_inq_format returned value",rv,"which is not recognized"))
}

