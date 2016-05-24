# Environment internal to the package
.bigalgebra_env <- new.env()

# Packgage library and global options
.onLoad <- function(libname, pkgname) {
  library.dynam("bigalgebra", pkgname, libname);
  options(bigalgebra.temp_pattern="matrix_")
  options(bigalgebra.mixed_arithmetic_returns_R_matrix=TRUE)
  options(bigalgebra.tempdir=tempdir)
  options(bigalgebra.DEBUG=FALSE)
}

.onUnload <- function(libpath) {
  library.dynam.unload("bigalgebra", libpath);
  options(bigalgebra.temp_pattern=c())
  options(bigalgebra.mixed_arithmetic_returns_R_matrix=c())
  options(bigalgebra.tempdir=c())
  options(bigalgebra.DEBUG=c())
}

