file_ext <- function(x) {
  # obtain the file extension 
  # copied from tools package 
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

obj.size.str <- function(x) {
  if (x >= 1024^3)       return(paste(round(x/1024^3, 1L), "Gb"))
  else if (x >= 1024^2)  return(paste(round(x/1024^2, 1L), "Mb"))
  else if (x >= 1024)    return(paste(round(x/1024, 1L), "Kb"))
  return(paste(x, "bytes")) 
} 

# not needed anymore since we have is.loaded in package:base,
# which is unknown to me before.
# 
#   is.null.ptr <- function(ns) {
#     .Call("is_Null_NS", ns) 
#   } 

## This is the implementation of 
## is_Null_NS in C, just in case we would need it in the future. 

##############################################################
#   #include <R.h>
#   #include <Rinternals.h>
#   
#   #ifdef __cplusplus
#   extern "C" {
#   #endif
#   
#   extern SEXP is_Null_NS(SEXP ns);
#   
#   #ifdef __cplusplus
#   }
#   #endif
#   
#   /*
#    * Tell if it is a NULL native symbol 
#    */
#   SEXP is_Null_NS(SEXP ns) {
#     SEXP ans;
#     PROTECT(ans = allocVector(LGLSXP, 1));
#     LOGICAL(ans)[0] = 1;
#     PROTECT(ns);
#     if (TYPEOF(ns) == EXTPTRSXP) { 
#       // Rprintf("ptr=%p.\n", EXTPTR_PTR(ns));
#       if (EXTPTR_PTR(ns) != NULL) LOGICAL(ans)[0] = 0;
#     } 
#     UNPROTECT(2);
#     return ans;
#   } 
#   
###############################################################

is.null.cxxfun <- function(cx) {
  # Tell if the returned object from cxxfunction in package inline
  # contains null pointer 
  env <- environment(cx@.Data)
  !is.loaded(env$f) 

  #  add <- body(cx@.Data)[[2]]
  #  # add is of class NativeSymbol
  #  .Call("is_Null_NS", add) 
} 

cxxfun.from.dll <- function(sig, code, DLL, check.dll = TRUE) { 
  # Create function objects from dll (most of the code are copied from
  # cxxfunction in package inline). 
  # 
  # Args:
  #  sig: a list of function signatures 
  #  DLL: object of class "DLLInfo"
  #  check.dll: check if the dll is loaded: When it is not 
  #    loaded, the function call might result in a segfault. 

  f <- DLL[['name']] 
  if (check.dll) {
    dlls <- getLoadedDLLs()
    if (!f %in% names(dlls)) 
      stop(paste("dso ", DLL[['path']], " is not loaded", sep = ''))
  } 

  res <- vector("list", length(sig))
  names(res) <- names(sig)
  res <- new("CFuncList", res)
 
  for(i in seq_along(sig)) {
    res[[i]] <- new("CFunc", code = code)
    fn <- function(arg) { NULL }

    ## Modify the function formals to give the right argument list
    args <- formals(fn)[rep(1, length(sig[[i]]))]
    names(args) <- names(sig[[i]])
    formals(fn) <- args

    ## create .Call function call that will be added to 'fn'
    body <- quote(.Call(EXTERNALNAME, ARG))[c(1:2, rep(3, length(sig[[i]])))]
    for (j in seq(along = sig[[i]])) body[[j + 2]] <- as.name(names(sig[[i]])[j])

    body[[1L]] <- .Call
    body[[2L]] <- getNativeSymbolInfo(names(sig)[[i]], DLL)$address
    ## update the body of 'fn'
    body(fn) <- body
    ## set fn as THE function in CFunc of res[[i]]
    res[[i]]@.Data <- fn
  }
  ## clear the environment
  rm(j)
  convention <- ".Call"
  if (identical(length(sig), 1L)) res[[1L]] else res
} 



cxxfun.from.dso.bin <- function(dso) {
  # Create function objects from dll (most of the code are copied from
  # cxxfunction in package inline). 
  # 
  # Args:
  #  dso: object of class cxxdso 
  # 
  # Note: we are assuming that the dso is not loaded so
  #   we create the dso file from the raw vector 
  #   and then loaded the dso. . 

  sig <- dso@sig 
  code <- dso@.MISC$cxxfun@code
  tfile <- tempfile() 
  f <- basename(tfile) 
  libLFile <- paste(tfile, ".", file_ext(dso@.MISC$dso.last.path), sep = '') 
  # write the raw vector containing the dso file to temporary file
  writeBin(dso@dso.bin, libLFile) 
  cleanup <- function(env) {
    if (f %in% names(getLoadedDLLs())) dyn.unload(libLFile)
      unlink(libLFile)
  }
  reg.finalizer(environment(), cleanup, onexit = TRUE)
  DLL <- dyn.load(libLFile) 
  assign('dso.last.path', libLFile, dso@.MISC) 
  res <- vector("list", length(sig))
  names(res) <- names(sig)
  res <- new("CFuncList", res)
  for(i in seq_along(sig)) {
    res[[i]] <- new("CFunc", code = code) 
    fn <- function(arg) { NULL }

    ## Modify the function formals to give the right argument list
    args <- formals(fn)[rep(1, length(sig[[i]]))]
    names(args) <- names(sig[[i]])
    formals(fn) <- args

    ## create .Call function call that will be added to 'fn'
    body <- quote(.Call(EXTERNALNAME, ARG))[c(1:2, rep(3, length(sig[[i]])))]
    for (j in seq(along = sig[[i]])) body[[j + 2]] <- as.name(names(sig[[i]])[j])

    body[[1L]] <- .Call
    body[[2L]] <- getNativeSymbolInfo(names(sig)[[i]], DLL)$address
    ## update the body of 'fn'
    body(fn) <- body
    ## set fn as THE function in CFunc of res[[i]]
    res[[i]]@.Data <- fn
  }
  ## clear the environment
  rm(j)
  rm(tfile) 
  convention <- ".Call"
  if (identical(length(sig), 1L)) res[[1L]] else res
} 


