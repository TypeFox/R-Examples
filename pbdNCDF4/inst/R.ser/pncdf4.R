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
  stop("--enable-parallel is no set. Please use nc_open or recompile pbdNCDF4.")
  # nc_open(filename, write = write, readunlim = readunlim, verbose = verbose)
} # End of nc_open_par().


nc_create_par <- function(filename, vars, force_v4 = TRUE, verbose = FALSE,
    comm = 0L, info = 0L){
  stop("--enable-parallel is no set. Please use nc_create or recompile pbdNCDF4.")
  # nc_create(filename, vars, force_v4 = force_v4, verbose = verbose)
} # End of nc_create_par().


nc_var_par_access <- function(nc, var, collective = TRUE, verbose = FALSE){
  # A serial dummy function.
  invisible()
} # End of nc_var_par_access().
