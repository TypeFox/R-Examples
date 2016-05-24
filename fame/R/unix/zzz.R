.onLoad <- function(libname, pkgname){
  ## this stuff needed to make Fame work
  ## just why Fame needs stuff from the Motif/LessTif library to work when
  ## there's no user interface involved is a question that will undoubtedly be
  ## pondered for the ages.
  libPath <- Sys.getenv("LD_LIBRARY_PATH")
  XmLib <- "/usr/lib/X11/../libXm.so"
  linkDir <- tempdir()
  tmpLink <- file.path(linkDir, "libXm.so.1")
  system(paste("ln -s", XmLib, tmpLink))
  Sys.setenv("LD_LIBRARY_PATH" = paste(libPath, linkDir, sep = ":"))
  Sys.setenv("FAME_EXPIRATION" = 0)
}
