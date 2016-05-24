nc_open_class_object <- function(filename, rv, write, readunlim, verbose){
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
} # End of nc_open_class_object().
