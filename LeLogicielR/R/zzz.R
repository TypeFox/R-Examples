.onLoad <- function(libname,pkg){

  library.dynam("LeLogicielR", pkg, libname)
  packageStartupMessage(utils::packageDescription('LeLogicielR',lib.loc = libname),appendLF=TRUE)

  y <- system.file("cor.test.2.sample.R",lib.loc = libname,package="LeLogicielR")
  source(y)
}
