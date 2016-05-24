#
# MigClim.genClust: Run the genetic clusters migration simulation.
#
MigClim.genClust <- function (hsMap="hsMap", barrier="barrier",
                              nrClusters=4, nrIterations=1, threshold=445,
                              outFile="out", initFile="")
{
  #
  # Get the number of rows and columns from the first input file.
  #
  Rst <- raster(paste(hsMap,"1.asc",sep=""))
  nrRows <- nrow(Rst)
  nrCols <- ncol(Rst)

  #
  # Call the genClust C function.
  #
  migrator <- .C("genClust", as.integer(nrRows), as.integer(nrCols),
                 as.integer(nrClusters), as.integer(nrIterations),
                 as.integer(threshold), hsMap, barrier, outFile,
                 initFile)
}

