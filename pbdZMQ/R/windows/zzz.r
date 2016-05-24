.onLoad <- function(libname, pkgname){
  dn <- paste(libname, "/", pkgname, "/libs",
              Sys.getenv("R_ARCH"), "/", sep = "")
  i.file <- "libzmq.dll"
  fn <- paste(dn, i.file, sep = "")
  ### Load "libzmq.dll".
  test <- try(dyn.load(fn, local = FALSE), silent = TRUE)
  if(class(test) == "try-error"){
    stop(paste("Could not load ", fn, sep = ""))
  }

  ### Load "pbdZMQ.dll".
  library.dynam("pbdZMQ", pkgname, libname)

  ### Preload to global environment.
  invisible(eval(parse(text = "pbdZMQ:::.zmqopt_init()")))

  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  ### Unload "pbdZMQ.dll".
  library.dynam.unload("pbdZMQ", libpath)

  dn <- paste(libpath, "/libs",
              Sys.getenv("R_ARCH"), "/", sep = "")
  i.file <- "libzmq.dll"
  fn <- paste(dn, i.file, sep = "")
  if(file.exists(fn)){
    ### Unload "libzmq.dll".
    test <- try(dyn.unload(fn), silent = TRUE)
  }

  invisible()
} # End of .onUnload().

.onAttach <- function(libname, pkgname){
  invisible()
} # End of .onAttach().

