#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Institut fuer Mathematik, Germany 2006
# *********************************************************
#
# Name:          readData.r
#                ---------------
# Author:        Joern Schulz
# Stand:         15.08.2006
#
#########################################################################
readData <- function(inputfile, DebugLevel="Normal")
{
# ========================================================================
#
# inputfile  (character)  'inputfile' giving the name of file to read, 
#            including pathname to the file. The extension of the name 
#            of file specified the format. Currently supported are ".txt", 
#            ".dat", ".dat.gz", ".fif", ".pet" and the graphics formats 
#            ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", 
#            ".jpg" and ".jpeg". See below to details to get more 
#            information about the formats.}
# DebugLevel (character)  This parameter controls the level of output. 
#            Following possibilities are available: The default "Normal" 
#            for standard level of output to screen or alternative "Detail" 
#            if it desirable to logged allmost all output to screen or 
#            "HardCore" for no information at all. }
#
# ========================================================================

 # =======================================================================
 # check the value 'inputfile'
  if(missing(inputfile)) {
      stop("The inputfile have to be of type character. It contains the path to filename and the filename. The path containing the path relatively to the working-directory of your R-session or containig the full path to the file.")
  }
  if(!is.character(inputfile)) 
      stop("'inputfile' have to be of type character.")
  if (length(inputfile) != 1)
      stop("Please specify exactly one 'inputfile'.")
  if (file.access(inputfile, 0))
      stop("The file, '", inputfile, "' doesn't exist.")
  if (file.access(inputfile, 4))
          stop("No read permission for file '", inputfile, "'.")

  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]  

 #######################################################
 # extraction of the file extension
  tmp <- unlist(strsplit(inputfile, ".", fixed=TRUE))
  fformat <- tmp[length(tmp)]
    

 # =======================================================================
 # main-routine
 #
  if (fformat == "txt"){
      fdata <- as.matrix(read.table(inputfile))

  } else if (fformat == "dat" || fformat=="gz"){
   ####################################################################
   #          same form of the header as pet but not binary

    fdata <- readASCIIData(inputfile, DebugLevel)
    
  } else if (fformat %in% c("fif", "pet")){
   ####################################################################
   #                             fif & pet
    
    fdata <- .Call("loadFile", inputfile, DebugLevel, PACKAGE="PET")

  } else if (fformat %in% c("tif", "tiff", "pgm", "ppm", "png", "pnm", "gif",
                            "jpg", "jpeg")){
   ####################################################################
   #     different graphic-formats supported by the adimpro library
    require(adimpro)
    fdata <- read.image(inputfile, compress=FALSE)$img
    #require(edges)
    #fdata <- read.image(inputfile)

  } else stop("The fformat ", fformat, " is not supported.")

  if (fformat == "fif")
      ans <- list(Signal=fdata[[1]], Header=fdata[2:10])
  else if (fformat == "pet")
      ans <- list(Signal=fdata[[1]], Header=fdata[2:5])
  else
      ans <- fdata


  return(ans)
}
  
