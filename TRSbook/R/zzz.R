.onLoad <- function(libname,pkg){
  library.dynam("TRSbook", pkg, libname)

  y <- system.file("cor.test.2.sample.R",lib.loc = libname,package="TRSbook")
  source(y)
}

