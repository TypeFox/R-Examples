#
# This .onLoad hook have been replaced by the '.registration' argument in
# useDynLib("Rmosek", .registration = TRUE)   // File: NAMESPACE
#
#.onLoad =
#function(libname, pkgname)
#{
#  dll <- library.dynam("Rmosek", package=pkgname, lib.loc=NULL) 
# 
#  syms = getNativeSymbolInfo(c("mosek_sym", "mosek_clean_sym", "mosek_version_sym", "mosek_read_sym", "mosek_write_sym"), dll)
#
#  # Create symbols in this package
#  env = environment(.onLoad)
#  sapply(names(syms), function(id) assign(id, syms[[id]], env))
#}


mosek = 
function(problem, opts=list())
{
  if (nargs() < 1) {
    print(help("mosek", package="Rmosek"))
    stop("Invalid number of arguments")
  }

  # Force evaluation of mandatory objects
  problem;

  r <- try(.Call(mosek_sym, problem, opts), silent = TRUE)

  if (inherits(r, "try-error")) {
    .Call(mosek_clean_sym)
    stop(r)
  }

  return(r)
}

mosek_clean = 
function()
{
  .Call(mosek_clean_sym) 
  return(invisible(NULL))
}

mosek_version =
function()
{
  r <- .Call(mosek_version_sym)
  return(r)
}

mosek_read = 
function(modelfile, opts=list())
{
  if (nargs() < 1) {
    print(help("mosek_read", package="Rmosek"))
    stop("Invalid number of arguments")
  }
 
  # Force evaluation of mandatory objects
  modelfile;
 
  r <- try(.Call(mosek_read_sym, modelfile, opts), silent = TRUE)

  if (inherits(r, "try-error")) {
    .Call(mosek_clean_sym)
    stop(r)
  }

  return(r)
}

mosek_write = 
function(problem, modelfile, opts=list())
{
  if (nargs() < 2) {
    print(help("mosek_write", package="Rmosek"))
    stop("Invalid number of arguments")
  }

  # Force evaluation of mandatory objects
  problem; modelfile;

  r <- try(.Call(mosek_write_sym, problem, modelfile, opts), silent = TRUE)

  if (inherits(r, "try-error")) {
    .Call(mosek_clean_sym)
    stop(r)
  }

  return(r)
}

