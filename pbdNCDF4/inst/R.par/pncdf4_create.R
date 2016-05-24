##GO## nc_create_preamble is an exact copy of top portion of nc_create
nc_create_preamble <- function(filename, vars, force_v4 = TRUE,
    verbose = FALSE){

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

        ###WCC: cut here.
        list(group = group, force_v4 = force_v4, vars = vars)
} # End of nc_create_preamble().


nc_create_finish <- function(nc, preamble, vars, verbose){
        ## restore preamble variables
        group <- preamble$goup

        ## nc_create_finish is an exact copy of nc_create section below this comment

        ### WCC: cut here.
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
} # End of nc_create_finish().
