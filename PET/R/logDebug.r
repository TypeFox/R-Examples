logDebug <- function(DebugLevel)
{
#
#
#

  DL1 <- FALSE
  DL2 <- FALSE
  if (is.character(DebugLevel) && DebugLevel=="Normal"){
      DL1 <- TRUE
  } else if (is.character(DebugLevel) && DebugLevel=="Detail"){
      DL1 <- TRUE
      DL2 <- TRUE
  } else if (is.character(DebugLevel) && DebugLevel=="HardCore"){
      # notthing is changed
  } else {
      cat("WARNING: DebugLevel='", DebugLevel,"' is not supported. \n", sep="")
      cat("Default 'Normal' is used. \n")
      DebugLevel <- "Normal"
      DL1 <- TRUE
  }

  invisible(list(DL1,DL2,DebugLevel))

}