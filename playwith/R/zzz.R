## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

## local environment to store state of plot windows
StateEnv <- new.env()

## local environment to store user options
.PlaywithEnv <- new.env()
.PlaywithEnv$options <- list()

.onLoad <- function(libname, pkgname)
{
    playwith.options(.defaultPlaywithOptions())
    ## set gWidgets toolkit (do not bother user)
    options(guiToolkit="RGtk2")
}
