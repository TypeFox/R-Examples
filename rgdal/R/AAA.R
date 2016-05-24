# following upgrade request 070122

if(!exists("Sys.setenv", envir = baseenv()))
     Sys.setenv <- Sys.putenv
.RGDAL_CACHE <- new.env(FALSE, parent=globalenv())
assign(".rgdal_old.PROJ_LIB", "", envir=.RGDAL_CACHE)
assign(".rgdal_old.GDAL_DATA", "", envir=.RGDAL_CACHE)

#.First.lib <- function(lib, pkg) {
.onLoad <- function(lib, pkg) {
#  require(methods, quietly = TRUE, warn.conflicts = FALSE)
#  require("sp")
  assign(".rgdal_old.PROJ_LIB", Sys.getenv("PROJ_LIB"), envir=.RGDAL_CACHE)
  assign(".rgdal_old.GDAL_DATA", Sys.getenv("GDAL_DATA"), envir=.RGDAL_CACHE)
  assign(".rgdal_old.NEEDED", FALSE, envir=.RGDAL_CACHE)
  if (file.exists(system.file("proj/nad.lst", package = "rgdal")[1])) {
    Sys.setenv("PROJ_LIB"=system.file("proj", package = "rgdal")[1])
    Sys.setenv("GDAL_DATA"=system.file("gdal", package = "rgdal")[1])
    assign(".rgdal_old.NEEDED", TRUE, envir=.RGDAL_CACHE)
  } else if (.Platform$OS.type == "windows") {
    assign(".rgdal_OSGeo4W", Sys.getenv("OSGEO4W_ROOT"), envir=.RGDAL_CACHE)
  }
  assign("OVERRIDE_PROJ_DATUM_WITH_TOWGS84", TRUE, envir=.RGDAL_CACHE)
  assign("silent", TRUE, envir=.RGDAL_CACHE)
  assign("has_proj_def.dat", .Call("PROJ4_proj_def_dat_Installed", PACKAGE="rgdal"), envir=.RGDAL_CACHE)

  library.dynam('rgdal', pkg, lib)

  .Call('RGDAL_Init', PACKAGE="rgdal")
}

.onAttach <- function(lib, pkg) {
  ver_ok <- getGDALCheckVersion()
  rver <- getGDALVersionInfo()

  gdl <- getGDAL_DATA_Path()
  pl <- getPROJ4libPath()
  if (nchar(pl) == 0) pl <- "(autodetected)"
  fn <- system.file("SVN_VERSION", package="rgdal")
  if (file.exists(fn)) {
    svn_version <- scan(fn, what=character(1), sep="\n", quiet=TRUE)
  } else {
    svn_version <- "(unknown)"
  }

  Smess <- paste('rgdal: version: ',
    utils::packageDescription("rgdal")$Version,
    ', (SVN revision ', svn_version, ')\n',
    ' Geospatial Data Abstraction Library ',
    'extensions to R successfully loaded\n',
    ' Loaded GDAL runtime: ', rver, ifelse(ver_ok, '\n',
    '\n   but rgdal build and GDAL runtime not in sync:\n   ... consider re-installing rgdal!!\n'),
    paste(" Path to GDAL shared files: ", gdl[1], sep=""), "\n",
    ifelse(GDAL_iconv(), "",
        paste(" GDAL does not use iconv for recoding strings.\n")),
    ' Loaded PROJ.4 runtime: ', getPROJ4VersionInfo(), '\n',
    paste(" Path to PROJ.4 shared files: ", pl[1], sep=""), "\n",
    ifelse(get("has_proj_def.dat", envir=.RGDAL_CACHE), "", "WARNING: no proj_defs.dat in PROJ.4 shared files\n"), sep="")
    splVersion <- version_sp_linkingTo()
  Smess <- paste(Smess, "Linking to sp version:", splVersion, "\n")
  spVcheck <- NULL
  if("sp" %in% .packages()) 
    spVcheck <- utils::packageVersion("sp") == splVersion
  if (!is.null(spVcheck) && !spVcheck) paste(Smess, 
    "sp version used to install rgdal and loaded sp version differ\n")
  packageStartupMessage(Smess, appendLF = FALSE)
}

#.Last.lib <- function(lib, pkg) {
.onUnload <- function(libpath) {
  if (get(".rgdal_old.NEEDED", envir=.RGDAL_CACHE)) {
    Sys.setenv("PROJ_LIB"=get(".rgdal_old.PROJ_LIB", envir=.RGDAL_CACHE))
    Sys.setenv("GDAL_DATA"=get(".rgdal_old.GDAL_DATA", envir=.RGDAL_CACHE))
  }
  .Call('RGDAL_Exit', PACKAGE="rgdal")
}

