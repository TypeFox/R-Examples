.onLoad <- function(libname, pkgname){
  ## If FAME is installed on this machine, load the real chli.dll first,
  ## then the fame.dll after.  Otherwise, load neither.

  ## first find out if we even have FAME
  fameDir <- Sys.getenv("FAME")
  if(nchar(fameDir) == 0){
    defaultDir <- "C:/Program Files/FAME"
    packageStartupMessage("FAME environment variable not set; assuming ", defaultDir)
    fameDir <- defaultDir
    fameDir <- sub("/$", "", chartr("\\", "/", fameDir))
  }
  
  if(!file.exists(fameDir)){
      packageStartupMessage("FAME directory ", fameDir, " not found.\n",
            "If you have FAME installed, specify it's location ",
            "via the FAME environment variable.\n",
            "Otherwise, this package is pretty useless.\n")
  }
  else {  ## apparently FAME is installed

    chliPath <- file.path(fameDir, "chli.dll")

    if(!file.exists(chliPath))
      stop(paste("chli.dll not found in", fameDir))
    dyn.load(chliPath, local = FALSE)
    if(is.loaded("cfmini"))
      assign(".chliPath", chliPath, pos = baseenv())
    else stop("Could not load chli.dll")
    
    library.dynam("fame", package = "fame", lib.loc = NULL)
    if(!is.loaded("dummyFameFunction"))
      stop("Could not load fame.dll")
  }
}

.onUnload <- function(libpath){
  if(exists(".chliPath", envir = baseenv())){
    dyn.unload(get(".chliPath", pos = baseenv()))
    remove(".chliPath", pos = baseenv())
    library.dynam.unload("fame", libpath)
  }
}

