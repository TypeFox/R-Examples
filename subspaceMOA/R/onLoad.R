.onLoad <- function(libname,pkgname) {
  #Reused code from the subspace package. This is done to keep WEKA from
  #creating a wekafiles folder in the user's home directory.
  if(is.na(isdir <-
             file.info(Sys.getenv("WEKA_HOME",
                                  file.path(path.expand("~"),
                                            "wekafiles")))$isdir) ||
       !isdir)
    Sys.setenv(WEKA_HOME = tempfile("RWeka"))
  
  rJava::.jpackage(name=pkgname,lib.loc=libname)
}