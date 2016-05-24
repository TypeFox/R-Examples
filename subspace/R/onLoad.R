.onLoad <- function(libname,pkgname) {
  #Load the jar files that this package uses.
  rJava::.jpackage(name=pkgname,lib.loc=libname)
  
  # Problem: loading the Weka jars creates a new directory called wekafiles in
  # the user's home directory. This is bad because nothing meaningful (for this
  # package) is saved in that directory and creating useless directories is Bad
  # Behavior and also against the policies of CRAN. 
  # Workaround: (proposed by Kurt Hornik)
  # If the directory given by
  # WEKA_HOME (or its default $HOME/wekafiles) was not created yet, make sure it
  # gets created in tempdir().
  if(is.na(isdir <-
             file.info(Sys.getenv("WEKA_HOME",
                                  file.path(path.expand("~"),
                                            "wekafiles")))$isdir) ||
       !isdir)
    Sys.setenv(WEKA_HOME = tempfile("subspace"))
}