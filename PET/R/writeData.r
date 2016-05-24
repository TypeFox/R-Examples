#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Institut fuer Mathematik, Germany 2006
# *********************************************************
#
# Name:          writeData.r
#                ---------------
# Author:        Joern Schulz
# Stand:         15.08.2006
#
#########################################################################
writeData <- function(data, outputfile, fileHeader=NULL, imType="normal",
                      fileOverwrite="ASK",DebugLevel="Normal")
{
# ========================================================================
#
# data     (matrix or array)  The data to be written out.
# outputfile (character)  'outputfile' naming the file to write to, including 
#          pathname to the file. The path containing the path relatively to the
#          working-directory of your R-session or containig the full path to the
#          file. The extension of the name of file specified the format. Currently
#          supported are ".txt", ".dat", ".dat.gz", ".fif", ".pet" and the graphics
#          formats ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" 
#          and ".jpeg". See below to details to get more information about the 
#          formats.
# fileHeader (list)  If file format ".dat", ".dat.gz", ".pet" or ".fif" then a 
#          fileheader is necessary. If in this cases \code{fileHeader=NULL} then
#          defaults values will be generate dependent on \code{imType}. See below 
#          to details to get more information.
# imType   (character)  'imType* specifies the type of image and with it the 
#          default values of \code{fileHeader}. Default is \code{imType="normal"}. 
#          Also implemented is \code{imType="radon"}, i.e. that the image is a 
#          sinogram (radon transformed image).
# fileOverwrite (character) Control the behaviour for overwriting a file. Supported
#          are "NO", "YES" and "ASK". If fileOverwrite = "ASK" the routine
#          will ask you in case of an existing outputfile with this name. Default
#          is "ASK".
# DebugLevel (character)  This parameter controls the level of output. Following
#          possibilities are available: The default "Normal" for standard level of
#          output to screen or alternative "Detail" if it desirable to logged 
#          allmost all output to screen or "HardCore" for no information at all.
#
# ========================================================================

 #######################################################
 # Check the value 'outputfile'
  if(missing(outputfile)) {
    stop("The outputfile have to be of type character. It contains the path to filename  and the filename. The path containing the path relatively to the working-directory  of your R-session or containig the full path to the file.")
  }
  if(!is.character(outputfile)) 
    stop("'outputfile' have to be of type character.")
  if (length(outputfile) != 1)
    stop("Please specify exactly one 'outputfile'.")
  if (is.na(match(fileOverwrite,c("ASK","YES","NO")))){
	 cat("WARNING: Parameter fileOverwrite= '", fileOverwrite, "' is not supported. Default parameter is used. \n",sep="")
	 fileOverwrite = "ASK";
  }
  if (!(file.access(outputfile, 0))){
	if (toupper(fileOverwrite) == "ASK"){
	  cat("WARNING: The file, '", outputfile, "' exist. If you like to \n",sep="")
	  cat("overwrite the file then type 'Yes' and enter. \n")
	  fileOverwrite <- readLines(n=1)
    } 
    if (!(toupper(fileOverwrite) == "YES")) 
        stop("The routine breaking off, beause it was given no permission to write.")
    if (file.access(outputfile, 2))
        stop("No write permission for file '", outputfile, "'.")
  }
  tmp <- unlist(strsplit(outputfile, "/", fixed=TRUE))
  if(length(tmp)==1){
    if (file.access(getwd(), 2))
        stop("No write permission in the working directory.")
  } else {
    tmp <- paste( tmp[-length(tmp)], collapse = "/" )
    if (file.access(tmp, 2))
        stop("No write permission in directory '", tmp, "'.")    
  }
  
  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]

 #######################################################
 # extraction of the file extension and of the filename without the pathname
  tmp <- unlist(strsplit(outputfile, ".", fixed=TRUE))
  fformat <- tmp[length(tmp)]
  fileName <- paste( tmp[-length(tmp)], collapse = "." )
  #tmp <- unlist(strsplit(tmp,"/", fixed=TRUE))
  #fileName <- tmp[length(tmp)]

 #######################################################
 #                   main-routine
  if (fformat == "txt"){
    fdata <- write.table(data,outputfile)

  } else if (fformat %in% c("dat","gz")){
   ####################################################################
   #           same form of the header as pet but not binary
    listNames <- names(fileHeader)
    header <- iniPetFifList(listType="pet", data=data, imType=imType)
    if("Description" %in% listNames)
        headerNames <- names(header)
    else {
        header <- list(SignalDim=header$SignalDim, 
                       XYmin=as.double(header$XYmin), 
                       DeltaXY=as.double(header$DeltaXY))
        headerNames <- names(header)
    }    
    corNames <- NULL
    incorNames <- NULL
    nonspecifyNames <- NULL
    corNames <- listNames[listNames %in% headerNames]
    incorNames <- listNames[!(listNames %in% headerNames)]
    nonspecifyNames <- headerNames[!(headerNames %in% corNames)]
    if (length(incorNames) != 0){
        if(DL1) cat("Following listnames are not supported: \n", incorNames,"\n")
    }
    if (length(nonspecifyNames) != 0){
        if(DL1) cat("Following listnames are not specified for a pet-file: \n", nonspecifyNames, "\n")
        if(DL1) cat("Default are used for the above names. \n")
    }
    header[corNames] <- fileHeader[corNames]
    
    if (fformat == "gz")
        writeASCIIData(data, fileName, fileHeader=header, Compress=TRUE, DebugLevel=DebugLevel)
    else
        writeASCIIData(data, outputfile, fileHeader=header, Compress=FALSE, DebugLevel=DebugLevel)
    
  } else if (fformat %in% c("fif", "pet", "mat")){
   ####################################################################
   #                     fif, mat & pet
    if (fformat == "fif" || fformat == "pet"){
       ################################################################
       #                        fif & pet
        listNames <- names(fileHeader)
        if (fformat == "pet")
            header <- iniPetFifList(listType="pet", data=data, imType=imType)
        else
            header <- iniPetFifList(listType="fif", data=data, imType=imType)
        headerNames <- names(header)
        corNames <- NULL
        incorNames <- NULL
        nonspecifyNames <- NULL
        corNames <- listNames[listNames %in% headerNames]
        incorNames <- listNames[!(listNames %in% headerNames)]
        nonspecifyNames <- headerNames[!(headerNames %in% corNames)]
        
        if (length(incorNames) != 0){
            if(DL1) cat("Following listnames are not supported: \n", incorNames,"\n")
        }
        if (length(nonspecifyNames) != 0){
            if (fformat == "pet"){
                if (DL1) cat("Following listnames are not specified for a pet-file: \n", nonspecifyNames, "\n")
            } else {
                if (DL1) cat("Following listnames are not specified for a fif-file: \n", nonspecifyNames, "\n")
            }
            if (DL1) cat("Default are uses for the above names. \n")
        }
        
        if (fformat == "fif" && ("FileName" %in% corNames) &&
                                fileHeader$FileName!=fileName)
            if (DL1) cat("WARNING: If fileHeader$FileName='", fileHeader$FileName, "' unequal to outputfile='", fileName, "' without file-typ, than default is use. \n", sep="")
        
        header[corNames] <- fileHeader[corNames]
        
        if (fformat == "pet"){
            CList <- list(SignalDim=header$SignalDim, 
                          Description=header$Description, 
                          XYmin=as.double(header$XYmin), 
                          DeltaXY=as.double(header$DeltaXY))
        } else {
            CList <- list(SignalDim=header$SignalDim, 
                          Description=header$Description, 
                          XYmin=as.double(header$XYmin), 
                          DeltaXY=as.double(header$DeltaXY),
                          SignalMinMax=as.double(header$SignalMinMax),
                          Date=header$Date,
                          FIFIdType=header$FIFIdType,
                          ArrayType=header$ArrayType)
        }
        
        .Call("writeFile", outputfile, as.double(data), CList, DebugLevel, PACKAGE="PET")
    
    } else {
       ###############################################################
       #                            mat
       CList <- list(SignalDim=dim(data))
       .Call("writeFile", outputfile, as.double(data), CList, DebugLevel, PACKAGE="PET")
    }

  } else if (fformat %in% c("tif", "tiff", "pgm", "ppm", "png", "pnm", "gif",
                              "jpg", "jpeg")){
   ####################################################################
   #     different graphic-formats supported by the adimpro library
    require(adimpro)
    data <- make.image(scaleImage(data, mode="max"), gamma=TRUE, compress=FALSE)
    write.image(data, file=outputfile)
    #require(edges)
    #write.image(img=data, file=outputfile, depth=16)

  } else stop("The fformat ", fformat, " is not supported.")

  invisible(NULL)
}
  
