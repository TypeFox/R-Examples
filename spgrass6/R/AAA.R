.GRASS_CACHE <- new.env(FALSE, parent=globalenv())
#assign("cmdCACHE", list(), envir=.GRASS_CACHE)
#assign("addEXE", .addexe(), envir=.GRASS_CACHE)

if(!exists("Sys.setenv", envir = baseenv())) Sys.setenv <- Sys.putenv

.onLoad <- function(lib, pkg) {
#.First.lib <- function(lib, pkg) {
#  if (Sys.getenv("OSTYPE") == "cygwin") {
#    cygwin_clean_temp(verbose=FALSE)
#    packageStartupMessage("Cygwin temp directory flushed\n", appendLF = FALSE)
#  }
# commented out for GRASS >= 6.3
#  require("sp")
  assign(".GRASS_old.GRASS_PAGER", Sys.getenv("GRASS_PAGER"), envir=.GRASS_CACHE)
  Sys.setenv("GRASS_PAGER"="cat")
  assign(".GRASS_old.GRASS_MESSAGE_FORMAT", Sys.getenv("GRASS_MESSAGE_FORMAT"),
    envir=.GRASS_CACHE)
  assign("INIT_USED", FALSE, envir=.GRASS_CACHE)
  
  Sys.setenv("GRASS_MESSAGE_FORMAT"="text")

  gisrc <- Sys.getenv("GISRC")
  loc <- Sys.getenv("LOCATION_NAME")

  assign("cmdCACHE", list(), envir=.GRASS_CACHE)
  assign("override_encoding", "", envir=.GRASS_CACHE)
  SYS <- ""
  if (.Platform$OS.type == "windows") {
    if (Sys.getenv("OSTYPE") == "msys") SYS <- "msys"
    else if (Sys.getenv("OSTYPE") == "cygwin") SYS <- "cygwin"
    else SYS <- "WinNat"
  } else if (.Platform$OS.type == "unix") SYS <- "unix"
  assign("SYS", SYS, envir=.GRASS_CACHE)
  res <- ""
  if (SYS == "msys" || SYS == "WinNat" || SYS == "cygwin") res =".exe"
  assign("addEXE", res, envir=.GRASS_CACHE)
  assign("WN_bat", "", envir=.GRASS_CACHE)
  if (SYS == "WinNat" && nchar(gisrc) > 0) {
    pyScripts <- basename(list.files(paste(Sys.getenv("WINGISBASE"),
      "scripts", sep="/"), pattern="py$"))
    names(pyScripts) <- sub("\\.py", "", pyScripts)
    assign("pyScripts", pyScripts, envir=.GRASS_CACHE)
  }

  assign("ignore.stderr", FALSE, envir=.GRASS_CACHE)
  assign("useGDAL", TRUE, envir=.GRASS_CACHE)
  assign("stop_on_no_flags_paras", TRUE, envir=.GRASS_CACHE)
  assign("plugin", NULL, envir=.GRASS_CACHE)
  assign("echoCmd", FALSE, envir=.GRASS_CACHE)
  assign("GV", "", envir=.GRASS_CACHE)
  assign("useIntern", FALSE, envir=.GRASS_CACHE)
  assign("legacyExec", .Platform$OS.type == "windows", envir=.GRASS_CACHE)
  assign("defaultFlags", NULL, envir=.GRASS_CACHE)
  assign("suppressEchoCmdInFunc", TRUE, envir=.GRASS_CACHE)
}

.onAttach <- function(lib, pkg) {
  gisrc <- Sys.getenv("GISRC")
  loc <- Sys.getenv("LOCATION_NAME")
  if (nchar(gisrc) == 0) gv <- "(GRASS not running)"
  else {
    gv <- .grassVersion() #system(paste("g.version", get("addEXE", envir=.GRASS_CACHE),sep=""),  intern=TRUE)
    comp <- .compatibleGRASSVersion(gv)
    if ( !comp ){
        stop( attr(comp, "message") )
    }
    assign("GV", gv, envir=.GRASS_CACHE)
    if(nchar(loc) == 0) {
      gisrc <- ifelse(.Platform$OS.type == "windows" &&
                (Sys.getenv("OSTYPE") == "cygwin"), 
		system(paste("cygpath -w", gisrc, sep=" "), intern=TRUE), 
		gisrc)
      loc <- read.dcf(gisrc)[1,"LOCATION_NAME"]
    }
  }

  Smess <- paste('GRASS GIS interface loaded ',
    'with GRASS version: ', gv, '\n',
    ifelse(nchar(loc) == 0, '', paste('and location: ', loc, '\n', sep="")),
      sep="")
  packageStartupMessage(Smess, appendLF = FALSE)
#  require("rgdal")
#  require("XML")
  
#  .GRASS_CACHE <- new.env(FALSE, parent=globalenv())
}

.onUnload <- function(lib, pkg) {
#.Last.lib <- function(lib, pkg) {
    if (get("INIT_USED", envir=.GRASS_CACHE)) {
        unlink_.gislock()
        unset.GIS_LOCK()
    }
    Sys.setenv("GRASS_PAGER"=get(".GRASS_old.GRASS_PAGER", envir=.GRASS_CACHE))
    Sys.setenv("GRASS_MESSAGE_FORMAT"=get(".GRASS_old.GRASS_MESSAGE_FORMAT",
        envir=.GRASS_CACHE))
    rm(.GRASS_CACHE)
}

