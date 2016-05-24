
library(pkgKitten)

setwd(tempdir())                        # which gets deleted anyway

kitten("simpleTest")                    # unfortunately testing fails :-/

R <- shQuote(file.path(R.home(component="bin"), "R"))

Sys.setenv("R_TESTS"="")                # needed for R CMD check; thanks for the tip, Hadley
Sys.setenv("R_LIBS"=paste(.libPaths(), collapse=.Platform$path.sep))

system(paste(R, "CMD build simpleTest"))

system(paste(R, "CMD check simpleTest_1.0.tar.gz"))



