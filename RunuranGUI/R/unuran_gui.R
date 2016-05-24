#############################################################################
##
##  Main
##
#############################################################################

## --------------------------------------------------------------------------
## First choice: use RGtk2

## options("guiToolkit"="RGtk2")

## --------------------------------------------------------------------------
## Main program

unuran.gui <- function(envir=parent.frame()) {

  ## check argument
  if (!is.environment(envir))
    stop ("Invalid argument 'envir'")
  
  ## main window
  main <- gwindow(title="Runuran")
  id(main) <- "Runuran"

  ## store environment where results are stored
  tag(main,"envir") <- envir

  ## Stage 1: Select type of distribution and generation method
  stage1(main)

  ## Stage 2 is started by button handler

  ## hide window
  invisible(main)
}

## --------------------------------------------------------------------------
