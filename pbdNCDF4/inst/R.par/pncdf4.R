#===============================================================================
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
##GO## added _par suffix and parameters comm, info
nc_open_par <- function(filename, write = FALSE, readunlim = TRUE,
    verbose = FALSE, comm = 0L, info = 0L){
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

        ##GO## added _par suffix - still need to add mpi parameters
	rv <- .C("R_nc4_open_par",
                 ## function parameters (copied, passed by address, returned as list)
                 as.character(filename),
                 as.integer(rv$cmode),
                 comm=as.integer(pbdMPI::comm.c2f(comm)),   ##GO##
                 info=as.integer(pbdMPI::info.c2f(info)),   ##GO##
                 id=as.integer(rv$id),	# note: nc$id is the simple integer ncid of the base file (root group in the file)
                 error=as.integer(rv$error),
                 ## .C parameters
                 PACKAGE="pbdNCDF4")
	if( rv$error != 0 )
		stop(paste("Error in nc_open trying to open file",filename))
	if( verbose )
		print(paste("nc_open: back from call to R_nc4_open, ncid=",rv$id))

        ###WCC: cut here.
        nc_open_class_object(filename, rv, write, readunlim, verbose)
} # End of nc_open_par().

nc_create_par <- function(filename, vars, force_v4 = TRUE, verbose = FALSE,
    comm = 0L, info = 0L){
        if(! force_v4){
          stop("This version will support parallel NetCDF version 4.")
        }

        ## call original nc_create preamble part
        preamble <- nc_create_preamble(filename, vars, force_v4 = force_v4,
                                       verbose = verbose )
        force_v4 <- preamble$force_v4
        vars <- preamble$vars

        ## middle section from nc_create with some modification
	nc <- list()

	#-----------------------------------------------------------
	# These values MUST MATCH the values in the source code file
	# ncdf.c, routine R_nc4_create!
	#-----------------------------------------------------------
	flag_NC_NOCLOBBER 	<- 1
	flag_NC_SHARE     	<- 2
	flag_NC_64BIT_OFFSET	<- 4
	flag_NC_NETCDF4		<- 8
        flag_NC_MPIIO           <- 16
        flag_NC_MPIPOSIX        <- 32

	#----------------
	# Create the file
	#----------------
	nc$cmode    <- 0
	if( force_v4 )
		nc$cmode <- nc$cmode + flag_NC_NETCDF4 + flag_NC_MPIIO
	nc$error    <- -1
	nc$id       <- -1
	if( verbose )
		print(paste("Calling R_nc4_create for file ",filename))
	nc<-.C("R_nc4_create_par",
               ## function parameters (copied, passed by address, retured as list)
		filename,
		as.integer(nc$cmode),
                comm=as.integer(pbdMPI::comm.c2f(comm)),  ##GO##
                info=as.integer(pbdMPI::info.c2f(info)),  ##GO##
		id=as.integer(nc$id),
		error=as.integer(nc$error),
               ## .C parameters
		PACKAGE="pbdNCDF4")
	if( nc$error != 0 )
		stop("Error in nc_create_par!")
	if( verbose )
		print(paste("back from R_nc4_create_par for file ",filename))
	nc$nvars  <- 0
	attr(nc,"class")  <- "ncdf4"
	nc$filename <- filename
	nc$writable <- TRUE

	nc$ndims  <- 0
	nc$dim    <- list()
	nc$var    <- list()

        return(nc_create_finish(nc, preamble, vars, verbose))
} # End of nc_create_par().


nc_var_par_access <- function(nc, var, collective = TRUE, verbose = FALSE){
  nc.format <- ncdf4_format(nc$id)
  if(nc.format != "NC_FORMAT_NETCDF4"){
    pbdMPI::comm.cat("\nnc.format = ", nc.format, " is not for parallel I/O.\n",
                     sep = "", quiet = TRUE)
    pbdMPI::comm.cat("Warning: Serial write need more caution in parallel.\n\n",
                     sep = "", quiet = TRUE)
    return(invisible())
  }

  if( verbose )
    print(paste("nc_var_par_access: entering with collective=",collective))

  if( class(nc) != "ncdf4" ) 
    stop("nc_var_par_access: passed nc NOT of class ncdf4!")

  if( verbose )
    print(paste("nc_var_par_access: ncid of file =",nc$id,
                "   filename=",nc$filename,"    writable=",nc$writable))

  if((mode(var) != 'character') &&
     (class(var) != 'ncvar4') &&
     (class(var) != 'ncdim4') &&
     (! is.na(var)))
    stop(paste("nc_var_par_access: 2nd argument must be an object of type ncvar",
               "(both parts of the ncdf object returned by nc_open()), the character-string name of a variable or dimension",
               "or NA to get the default variable from the file.  If the file is netcdf version 4",
               "format and uses groups, then the fully qualified var name must be given, for",
               "example, model1/run5/Temperature"))

  varid <- vobjtovarid4( nc, var, verbose=verbose, allowdimvar=TRUE )$id


  ncid <- nc$id

  if(collective)
    access <- 1
  else
    access <- 0

  error <- -1

  if(verbose){
      print(c("str(ncid)", "str(varid)", "str(access)", "str(error)"), "all")
  }
  
  nc<-.C("R_nc4_var_par_access",
         ## function parameters (copied, passed by address, returned as list)
         as.integer(ncid),
         as.integer(varid),
         as.integer(access),
         error=as.integer(error),
         ## .C parameters
         PACKAGE="pbdNCDF4")

  if( nc$error != 0 )
    stop("Error in nc_var_par_access!")
  if( verbose )
    print("back from R_nc4_var_par_access")
  invisible(nc$error)
} # End of nc_var_par_access().
