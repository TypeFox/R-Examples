### Lastest load into a package.

### Export Namespace does not use .First.lib() and .Last.lib(), but use
### .onLoad() and .onUnload().
# .First.lib <- function(lib, pkg){
# } # End of .First.lib().

# .Last.lib <- function(libpath){
# } # End of .Last.lib().

.onLoad <- function(libname, pkgname){
  library.dynam("pmclust", pkgname, libname)
  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  library.dynam.unload("pmclust", libpath)
  invisible()
} # End of .onUnload().
