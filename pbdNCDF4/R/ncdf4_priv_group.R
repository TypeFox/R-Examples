
#===============================================================
# Given the integer root_id, returns the integer IDs of the
# groups BELOW this node.  It is entirely possible that the
# return value will be an array of length 0, if there are no
# nodes below this one.
#
nc_grpids <- function( root_id ) {

	if( ! is.integer(root_id))
		root_id <- as.integer( root_id + 0.1 )

	rv <- list()
	rv$root_id <- root_id
	rv$ngrps   <- -1
	rv$error   <- -1
	rv <- .C("R_nc4_inq_ngroups",
		as.integer(rv$root_id),
		ngrps=as.integer(rv$ngrps),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop(paste("Error in nc_grpids trying to get number of groups; C function R_nc4_inq_ngroups returned erro"))

	ngrps <- rv$ngrps
	if( ngrps == 0 )
		return		# if NO groups, return nothing

	#-----------------------------------------------------
	# Now that we have the number of groups, get their IDs
	#-----------------------------------------------------
	rv         <- list()
	rv$root_id <- root_id
	rv$gids    <- integer(ngrps)
	rv$error   <- -1
	rv <- .C("R_nc4_inq_groupids", 
		as.integer(rv$root_id),
		gids=as.integer(rv$gids),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop("C function R_nc4_inq_groupids returned error")

	return( rv$gids )
}

#===============================================================
# Given the ID of a group, returns that group's name
#
nc_grpname <- function( root_id ) {

	if( ! is.integer(root_id))
		root_id <- as.integer( root_id + 0.1 )

	ierr <- as.integer(-1)
	rv <- .Call( "R_nc4_grpname", as.integer(root_id), ierr, PACKAGE="pbdNCDF4" )
	if( ierr != 0 ) {
		stop(paste("Error in nc_grpname, encountered when root_id=",root_id))
		}
	return( rv )
}

#===================================================================================================
# Fill out information for a group, given its id (which can be the ID of the root group, 
# of course!).  Note that the fqgn of the ROOT group is just the empty string, "".
# 
# This sets:
#	$name: 		group name, or "" for the root group. Example: "run1"
#	$fqgn: 		fully qualified group name, or "" for the root group. Example: "model1/run1"
#	$id:		the group id (equal to the ncid for the root group)
#	$nvars:		number of variables in this group
#	$ndims:		number of dimensions in this group (note: vars can use dimensions
#			in this group OR ANY PARENT GROUP)
#	$dimid:		Dimids of the dims that are visible to this group
#	$natts:		number of attributes in this group
#
nc_get_grp_info <- function( gid, parent_fqgn, format ) {	# sets $name, $fqgn, $id, $nvars, $ndims, $natts

	retval <- list( id=gid )

	if( format == 'NC_FORMAT_NETCDF4' )
		retval$name  <- nc_grpname( gid )
	else
		retval$name = ""

	#----------------------------------------
	# Get standard stuff: ndims, nvars, natts
	#----------------------------------------
	rv       <- list()
	rv$ndims <- -1
	rv$nvars <- -1
	rv$natts <- -1
	rv$error <- -1
	rv <- .C("R_nc4_inq",
		as.integer(gid),
		ndims=as.integer(rv$ndims),
		nvars=as.integer(rv$nvars),
		natts=as.integer(rv$natts),
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop(paste("nc_get_grp_info: R_nc4_inq returned error on group id", gid ))

	retval$ndims <- rv$ndims
	retval$nvars <- rv$nvars
	retval$natts <- rv$natts

	#----------------------------------------
	# Get the dim IDs that this group uses.
	# In v4, it's no longer simply 0..ndims-1
	#----------------------------------------
	rv              <- list()
	rv$group_dimids <- array( as.integer(-1), retval$ndims )
	rv$error        <- -1
	rv <- .C("R_nc4_inq_dimids",
		as.integer(gid),
		group_dimids=rv$group_dimids,
		error=as.integer(rv$error),
		PACKAGE="pbdNCDF4")
	if( rv$error != 0 ) 
		stop(paste("nc_get_grp_info: R_nc4_inq_dimids returned error on group id", gid ))

	retval$dimid <- rv$group_dimids
	if( sum( is.na(retval$dimid) ) > 0 ) {
		print(paste("Error getting group info for group", retval$name,":"))
		print(paste("got some NA's as the dimids.  There are ", retval$ndims, " dims in this group,"))
		print("and here are the IDs:")
		print(retval$dimid)
		stop("cannot have NA's as dimids!")
		}

	#--------------------------------------------------------
	# Get the fully qualified group name (fqgn).  This is ""
	# for the root group, something like "model1" for a
	# group one level down, for 2 levels "model1/run1", etc.
	#--------------------------------------------------------
	if( retval$name == "/" ) {
		retval$name  <- ""
		retval$fqgn  <- ""
		}
	else
		{
		if( parent_fqgn == "" )
			retval$fqgn  <- retval$name
		else
			retval$fqgn  <- paste( parent_fqgn, "/", retval$name, sep='' )
		}

	class( retval ) <- "ncgroup4"

	return( retval )
}

#===============================================================
# Given a "group" (list with components $name, $fqgn, and $id), this returns the group list
# of all the groups BELOW this group.
# Return val is a list (possibly of length 0!) with entries $name, $fqgn, $id, $nvars, $ndims, $natts
#
nc_groups_below <- function( root_group, format ) {

	if( class(root_group) != "ncgroup4" )
		stop(paste("error, nc_groups_below must be called with a ncgroup-class object"))

	root_id <- root_group$id

	retval <- list()
	nret   <- 0

	#-----------------------------------
	# Get IDs of groups below this group
	#-----------------------------------
	gids <- nc_grpids( root_id ) 	# returns array of integers (possibly length 0!) of IDs of groups BELOW this group
	n_below <- length( gids )

	if( n_below == 0 )
		return( retval )	# returns an empty list

	for( ib in 1:n_below ) {
		newgroup <- nc_get_grp_info( gids[ib], root_group$fqgn, format )
		nret  <- nret + 1
		retval[[nret]] <- newgroup

		gg <- nc_groups_below( retval[[nret]], format )
		nn <- length(gg)
		for( ii in nc4_loop(1,nn) ) {
			nret <- nret + 1
			retval[[nret]] <- gg[[ii]]
			}
		}

	return( retval )
}

#===============================================================
# This returns a list of unique groups.  Each group has the following
# elements:  
# 	$name		(group name.  Has NO slashes. Note: may be duplicated, since group names are not necessarily unique!)
# 	$fqgn		(fully qualified group name. Does NOT start or end with a slash. This is unique across a file. "names(list)" is also set to this, for convenience)
# 	$fqpn		(fully qualified parent group name.  Does NOT end with a slash.)
#	$var[[]] 	(list of ncvar objects) 
#	$dim[[]] 	(list of ncdim objects)
#
# The alg guarantees that the order of the groups in the returned
# list goes from the root outward; i.e., by the time some group
# is encountered on the list, it is guaranteed that all of that
# group's ancestors will already have been encountered on the list.
#
nc_parse_group_structure <- function( vars, verbose=FALSE ) {

	debug <- verbose
	if( debug ) print('nc_parse_group_structure: entering')

	group <- list()

	root_grp <- list( name='', fqgn='', fqpn='', var=list(), dim=list() )  # Note: fqgn is a fully qualified group name, i.e., /group1/model1 etc.
	group[[1]] <- root_grp
	names(group)[1] <- ''	# set 'names' array to fqgn

	#----------------------------------------------------------------------------------
	# Go through each variable.  See if it lives in a specified group.  If it does NOT, 
	# then add this var to the root group.  If it DOES, then make sure that group is
	# on the group list, and add this var to that group.
	#----------------------------------------------------------------------------------
	nv <- length(vars)
	for(iv in 1:nv) {

		vv <- vars[[iv]]

		if( debug ) print(paste('------------- nc_parse_group_structure: working on var', vv$name))

		if( nslashes_ncdf4( vv$name ) == 0 ) {		# look at FULLY QUALIFIED var name, which should have exacty zero leading forward slashes if at root group
			#---------------------------------------------------
			# This var has no group, so put it in the root group
			#---------------------------------------------------
			if( debug ) print(paste('nc_parse_group_structure: var', vv$name, 'is in root group'))
			#----------------------------------------------------------------------
			# Make sure there is not already a var with this name in the root group
			#----------------------------------------------------------------------
			nv <- length(group[[1]]$var)
			for( iv in nc4_loop(1,nv) ) {
				if( vv$name == group[[1]]$var[[iv]]$name ) {
					stop(paste("Error, trying to add var named", 
						vv$name, " to the root group, but there is already ",
						"a var named ", group[[1]]$var[[iv]]$name,
						"in the root group.  Have you mistakenly added ",
						"the same var twice?  Have you correctly specifed ",
						"the groups of any variables with duplicated names?",
						"Syntax for indicating group is a var named, for example, /groupname/varname"))
					}
				}
			group[[1]]$var[[nv+1]] <- vv
			}
		else
			{
			#--------------------------------------------------------------------
			# This var has a group, see if its group is already on the group list
			#--------------------------------------------------------------------
			ng <- length(group)
			my_fqgn <- nc4_basename( vv$name, dir=TRUE )
			if( debug ) print(paste('nc_parse_group_structure: var', vv$name, 'is in group', my_fqgn))
			gidx <- -1
			for( ig in 1:ng ) {
				if( my_fqgn == group[[ig]]$fqgn ) {
					gidx <- ig
					break
					}
				}
			if( gidx == -1 ) {
				#---------------------------------------------
				# This group is not on the list yet, so add it
				#---------------------------------------------
				if( debug ) print(paste('nc_parse_group_structure: adding group', my_fqgn, ' to group list as a new entry'))
				ng <- length(group)
				my_fqpn <- nc4_basename(my_fqgn,dir=TRUE)
				group[[ng+1]] <- list( name=nc4_basename(my_fqgn), fqgn=my_fqgn, fqpn=my_fqpn, var=list( vv ), dim=list() )
				names(group)[ng+1] <- my_fqgn		# set 'names' array to fqgn
				}
			else
				{
				#-------------------------------------------------
				# This group is already on the list at index gidx; 
				# add this var to that list
				#-------------------------------------------------
				if( debug ) print(paste('nc_parse_group_structure: group', my_fqgn, ' is on the group list, so just adding var', vv$name,' to that group on the list'))
				#-----------------------------------------------------------------
				# Make sure there is not already a var with this name in the group
				#-----------------------------------------------------------------
				nv <- length( group[[gidx]]$var )
				for( iv in nc4_loop(1,nv) ) {
					if( vv$name == group[[gidx]]$var[[iv]]$name ) {
						stop(paste("Error, trying to add var named", 
							vv$name, " to group", group[[gidx]]$name, "but there is already ",
							"a var named ", group[[gidx]]$var[[iv]]$name,
							"in that group.  Have you mistakenly added ",
							"the same var twice?  Have you correctly specifed ",
							"the groups of any variables with duplicated names?"))
						}
					}
				group[[gidx]]$var[[nv+1]] <- vv
				}
			}

		#------------------------------------------------------------
		# Now go through this var's dims and see where they will live
		#------------------------------------------------------------
		nd <- vv$ndims
		for( idim in nc4_loop(1,nd) ) {
			dd <- vv$dim[[idim]]
			if(debug) print(paste("working on var",vv$name,"'s dim",idim,"named",dd$name))
			if( nslashes_ncdf4( dd$name ) == 0 ) {		# look at FULLY QUALIFIED dim name, which should have exacty zero leading forward slashes if at root group
				#---------------------------------------------------
				# This dim has no group, so put it in the root group
				#---------------------------------------------------
				if( debug ) print(paste('nc_parse_group_structure: dim', dd$name, 'has no specified group (and so is in root group)'))
				#-------------------------------------------------------
				# If there is already a dim with this name in the group,
				# make sure it's the SAME dim.
				#-------------------------------------------------------
				nd <- length(group[[1]]$dim)
				found_dim_already <- FALSE
				for( id in nc4_loop(1,nd) ) {
					if( dd$name == group[[1]]$dim[[id]]$name) {
						#-----------------------------------------------
						# HAVE seen this dim; ensure it's the same dim!!
						#-----------------------------------------------
						if(! ncdim_same( dd, group[[1]]$dim[[id]])) {
							stop(paste("Error, trying to add dim named", 
								dd$name, " to the root group, but there is already ",
								"a dim named ", group[[1]]$dim[[id]]$name,
								"in the root group, with different length.  Have you mistakenly resued ",
								"a dimension name for two different dims?  Have you correctly specifed ",
								"the groups of any dims with duplicated names?"))
							}
						found_dim_already <- TRUE
						if( debug ) print(paste("dim",dd$name,"is already in root group dim list, entry number", id))
						}
					}
				if( ! found_dim_already ) {
					#-----------------------------------------------
					# Only add this dim if it's not been seen before
					#-----------------------------------------------
					if( debug ) print(paste("adding dim",dd$name,"to root group dim list, entry number", nd+1))
					group[[1]]$dim[[nd+1]] <- dd
					}
				}
			else
				{
				#--------------------------------------------------------------------
				# This dim has a group, see if its group is already on the group list
				#--------------------------------------------------------------------
				ng <- length(group)
				my_fqgn <- nc4_basename( dd$name, dir=TRUE )
				if( debug ) print(paste('nc_parse_group_structure: dim', dd$name, 'is in fully qualified group', my_fqgn))
				gidx <- -1
				for( ig in 1:ng ) {
					if( my_fqgn == group[[ig]]$fqgn ) {
						gidx <- ig
						break
						}
					}
				if( gidx == -1 ) {
					#---------------------------------------------
					# This group is not on the list yet, so add it
					#---------------------------------------------
					if( debug ) print(paste('nc_parse_group_structure: adding (dim) group', my_fqgn, ' to group list as a new entry'))
					ng <- length(group)
					my_fqpn <- nc4_basename(my_fqgn,dir=TRUE)
					group[[ng+1]] <- list( name=nc4_basename(my_fqgn), fqgn=my_fqgn, fqpn=my_fqpn, var=list(), dim=list(dd) )
					names(group)[ng+1] <- my_fqgn		# set 'names' array to fqgn
					}
				else
					{
					#-------------------------------------------------
					# This group is already on the list at index gidx; 
					# add this dim to that list
					#-------------------------------------------------
					if( debug ) print(paste('nc_parse_group_structure: (dim"s) group', my_fqgn, ' is on the group list, so just adding dim', dd$name,' to that group on the list'))
					#---------------------------------------------------------------------------
					# Make sure there is not already a different dim with this name in the group
					#---------------------------------------------------------------------------
					nd <- length( group[[gidx]]$dim )
					found_dim_already <- FALSE
					for( id in nc4_loop(1,nd) ) {
						if( dd$name == group[[gidx]]$dim[[id]]$name ) {
							if( ! ncdim_same( dd, group[[gidx]]$dim[[id]] ))
								stop(paste("Error, trying to add dim named", 
									dd$name, " to group", group[[gidx]]$name, "but there is already ",
									"a dim named ", group[[gidx]]$dim[[id]]$name,
									"in that group.  Have you mistakenly reused ",
									"the same dim name with different dims?  Have you correctly specifed ",
									"the groups of any dims with duplicated names?"))
							found_dim_already <- TRUE
							}
						}
					if( ! found_dim_already ) {
						if( debug ) print(paste("adding dim",dd$name,"to group", 
							group[[gidx]]$name, "dim list, entry number", nd+1))
						group[[gidx]]$dim[[nd+1]] <- dd
						}
					}
				}
			}
		}

	#-------------------------------------------------------------------------
	# Now calculate the level (directory depth) of each group. The root group
	# is always at level 1.  Any groups under the root group are at level
	# 2.  Groups under that are level 3.  Etc.
	# Note that this code assumes:
	#	1) root group is named ""
	#	2) the names of groups under the root group neither start nor
	#	   end with a slash.  Examples: "model1", "model1/run1", etc.
	#-------------------------------------------------------------------------
	ng <- length( group )
	maxlevel <- -1
	for( ig in 1:ng ) {
		nam <- group[[ig]]$fqgn 
		if( nam == "" )		# root group
			group[[ig]]$level <- 1
		else
			group[[ig]]$level <- nslashes_ncdf4( nam ) + 1
		if( group[[ig]]$level > maxlevel ) maxlevel <- group[[ig]]$level
		}

	#-----------------------------------------------------------------------
	# Rearrange the group list so that levels go from least to greatest
	#-----------------------------------------------------------------------
	tlist <- list()
	idx_out <- 0
	for( ilev in 1:maxlevel ) {
		for( ii in 1:ng ) {
			if( group[[ii]]$level == ilev ) {
				idx_out <- idx_out + 1
				tlist[[idx_out]] <- group[[ii]]
				}
			}
		}
	group <- tlist

	if( debug ) {
		print("----------- group list: ----------" )
		ng <- length(group)
		for( ig in 1:ng ) {
			print(paste("Group ", ig, ", named ", group[[ig]]$fqgn,
				", is at level ", group[[ig]]$level,
				" and has full qual parent named ", group[[ig]]$fqpn,
				":", sep='' ))
			nd <- length( group[[ig]]$dim )
			for( id in nc4_loop(1,nd) ) 
				print(paste( ".... dim ", id, ": ", group[[ig]]$dim[[id]]$name, 
					' (length=', group[[ig]]$dim[[id]]$len, ')', sep='' ))
			nv <- length( group[[ig]]$var )
			for( iv in nc4_loop(1,nv) ) 
				print(paste( ".... var ", iv, ": ", group[[ig]]$var[[iv]]$name, sep='' ))
			}
		}

	return( group )
}

#===================================================================
# Inputs:
#	parentid: simple integer ID of parent of the group to be created
#	group_name: name of the group to be created
#
# Return value:
#	simple integer group ID of the new group
#
# Routine halts on error
#
nc_make_group_inner <- function( parentid, group_name ) {

	if( ! is.numeric(parentid))
		stop(paste("Error, first input (parentid) must be a simple integer"))

	if( ! is.character(group_name))
		stop(paste("Error, second input (group_name) must be a character string"))

	newgroup <- list( gid=-1, error=-1 )

	newgroup <-.C("R_nc4_def_grp",
		as.integer(parentid),
		as.character(group_name),
		gid=as.integer(newgroup$gid),
		error=as.integer(newgroup$error),
		PACKAGE="pbdNCDF4")

	if( newgroup$error != 0 )
		stop(paste("Error in R_nc4_def_grp trying to make group named", group_name))

	return( newgroup$gid )
}

