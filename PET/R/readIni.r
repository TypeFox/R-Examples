readIni <- function(fileName, DebugLevel="Normal")
{
#
#
#

###############################################################
#          check the value 'fileName' and DebugLevel
if(missing(fileName)) {
    stop("The fileName has to be of type character. It contains the path to filename and the filename. The path containing the path relatively to the working-directory of your R-session or containig the full path to the file.")
    }
if(!is.character(fileName)) 
    stop("'fileName' have to be of type character.")
if (length(fileName) != 1)
    stop("Please specify exactly one 'fileName'.")
if (file.access(fileName, 0))
    stop("The file, '", fileName, "' doesn't exist.")
if (file.access(fileName, 4))
    stop("No read permission for file '", fileName, "'.")

DL1 <- logDebug(DebugLevel)
DebugLevel <- DL1[[3]]
DL2 <- DL1[[2]]
DL1 <- DL1[[1]]

###############################################################
#                      Open connection
if (DL1)
    cat("Open file: ",fileName, "\n")

con <- file(fileName, "rb")
on.exit(close(con)) # be careful ...

###############################################################
#               Reading the INI-File
if (DL2)
    cat("Reading INI-file.\n")

tmp <- readLines(con, n=1)
tmp <- unlist(strsplit(tmp, "=", fixed=TRUE))
tmp <- c(tmp[1],paste(tmp[2:length(tmp)], collapse="=" ))
if (!(tmp[1]=="mode")){
    cat("The first line of the INI-file must contain the parameter mode and not '", tmp[1],"'. \n", sep="")
    cat("See help to get more information about the structure of the file. \n")
    stop()
} else {
    mode <- tmp[2]
    if (DL2)
        cat("'mode' entry is set to '", mode, "' \n", sep="")
}


if (mode %in% c("FB", "BF", "CNF", "CBF", "CC", "CNC", "CBC", "Test", "ART", "EM", "CG")){

if (mode %in% c("FB", "BF", "CNF", "CBF", "CC", "CNC", "CBC", "Test")){
    ###########################################################
    #        direct reconstruction-method are used

    # Initialisation a list with default values
    leIniList <- 12
    
    iniList <- list()
    length(iniList) <- leIniList
    names(iniList) <- c("mode", "rData", "XSamples", "YSamples", "Xmin", "Ymin", 
           "DeltaX", "DeltaY", "InterPol", "FilterTyp", "DebugLevel", "oData")
    
    iniList$mode <- mode
    iniList$InterPol <- 1
    iniList$FilterTyp <- "Hamming"
    iniList$DebugLevel <- DebugLevel
        
    listNames <- names(iniList)
    listNamesChar <- c("rData", "FilterTyp", "DebugLevel", "oData")
    listNamesNumeric <- listNames[!(listNames %in% listNamesChar)]

} else {
    ###########################################################
    #        iterative reconstruction-method are used
    
    # Initialisation a list with default values
    leIniList <- 32
    
    iniList <- list()
    length(iniList) <- leIniList
    names(iniList) <- c("mode", "rData", "XSamples", "YSamples", "StartImage", "UseFast",
          "RadonKernel", "Iterations", "IterationsType", "SaveIterations",
          "SaveIterationsName", "LowestALevel", "ConstrainMin", "ConstrainMax",
          "Alpha", "Beta", "Regularization", "KernelFileSave", 
          "KernelFileName", "RefFileName", "ThetaSamples", "RhoSamples",
          "ThetaMin", "RhoMin", "DeltaTheta", "DeltaRho", "Xmin", "Ymin", 
          "DeltaX", "DeltaY", "OverSamp", "DebugLevel")
    
    iniList$mode <- mode
    iniList$StartImage <- "None"
    iniList$UseFast <- 1
    iniList$RadonKernel <- "NN"
    iniList$Iterations <- 20
    iniList$IterationsType <- "random"
    iniList$SaveIterations <- 0
    iniList$SaveIterationsName <- ""
    iniList$LowestALevel <- 0.0
    iniList$ConstrainMin <- -1
    iniList$ConstrainMax <- -1
    iniList$Alpha <- 1.0
    iniList$Beta <- 1.0
    iniList$Regularization <- 0
    iniList$KernelFileSave <- 0
    iniList$KernelFileName <- ""
    iniList$RefFileName <- "None"
    iniList$OverSamp <- 0
    iniList$DebugLevel <- DebugLevel 
    
    listNames <- names(iniList)
    listNamesChar <- c("rData", "StartImage", "RadonKernel", "IterationsType",
                       "SaveIterationsName", "KernelFileName", "RefFileName",
                       "DebugLevel")
    listNamesNumeric <- listNames[!(listNames %in% listNamesChar)]
}
    
tmp <- readLines(con, n=leIniList)
    
for (i in 1:length(tmp)){
    tmp2 <- unlist(strsplit(tmp[i], "=", fixed=TRUE))
    tmp2 <- c(tmp2[1],paste(tmp2[2:length(tmp2)], collapse="=" ))
    if (tmp2[1] %in% listNames){
        if (tmp2[1] %in% listNamesNumeric){
            iniList[tmp2[1]] <- as.numeric(tmp2[2])
            if (DL2) cat(tmp2[1], "entry is set to", tmp2[2], "\n")
        } else {
            iniList[tmp2[1]] <- tmp2[2]
            if (DL2) cat(tmp2[1], " entry is set to '", tmp2[2], "'\n", sep="")
        }
    } else {
        cat("WARNING: '", tmp2[1], "' is an illegal parameter and is ignored. \n", sep="")
    }
}


if (mode %in% c("FB", "BF", "CNF", "CBF", "CC", "CNC", "CBC", "Test")){
    ###########################################################
    #  parameter checking for iterative reconstruction-method
    
    # load rData
    if (is.null(iniList$rData))
         stop("'rData' entry missing in INI file.")
    tmp <- unlist(strsplit(iniList$rData, ".", fixed=TRUE))
    fformatRadon <- tmp[length(tmp)]
    rfile <- iniList$rData
    if (DL2) cat("\n")
    iniList$rData <- readData(iniList$rData, DebugLevel)$Signal
    if (DL2) cat("\n")
    
    # If oData specified than loading oData
    if (!(is.null(iniList$oData))){
        tmp <- unlist(strsplit(iniList$oData, ".", fixed=TRUE))
        fformatOrginal <- tmp[length(tmp)]
        sfile <- iniList$oData
        iniList$oData <- readData(iniList$oData, DebugLevel)
        if (DL2) cat("\n")
    }
    
    ##############################################
    # checking the following parameters:
    # XSamples, YSamples, XMin, YMin, DeltaX and DeltaY
    if ( !(is.null(iniList$oData)) && fformatOrginal %in% c("pet", "fif", "dat")){
        if(DL1) cat("Setting unknown parameters from '", sfile, "' \n", sep="")
        # check XSamples
        if(is.null(iniList$XSamples)){
            iniList$XSamples <- iniList$oData$Header$SignalDim[[1]]
            if(DL2)
                cat("'XSamples' entry is set to", iniList$XSamples, "\n")
        } else if (iniList$oData$Header$SignalDim[[1]] != iniList$XSamples){
            if(DL2)
                cat("Note: 'XSamples' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check YSamples
        if(is.null(iniList$YSamples)){
            iniList$YSamples <- iniList$oData$Header$SignalDim[[2]]
            if(DL2)
                cat("'YSamples' entry is set to", iniList$YSamples, "\n")
        } else if (iniList$oData$Header$SignalDim[[2]] != iniList$YSamples){
            if(DL2)
                cat("Note: 'YSamples' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check Xmin
        if(is.null(iniList$Xmin)){
            iniList$Xmin <- iniList$oData$Header$XYmin[[1]]
            if(DL2)
                cat("'Xmin' entry is set to", iniList$Xmin, "\n")
        } else if (iniList$oData$Header$XYmin[[1]] != iniList$Xmin){
            if(DL2)
                cat("Note: 'Xmin' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check Ymin
        if(is.null(iniList$Ymin)){
            iniList$Ymin <- iniList$oData$Header$XYmin[[2]]
            if(DL2)
                cat("'Ymin' entry is set to", iniList$Ymin, "\n")
        } else if (iniList$oData$Header$XYmin[[2]] != iniList$Ymin){
            if(DL2)
                cat("Note: 'Ymin' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check DeltaX
        if(is.null(iniList$DeltaX)){
            iniList$DeltaX <- iniList$oData$Header$DeltaXY[[1]]
            if(DL2)
                cat("'DeltaX' entry is set to", iniList$DeltaX, "\n")
        } else if (iniList$oData$Header$DeltaXY[[1]] != iniList$DeltaX){
            if(DL2)
                cat("Note: 'DeltaX' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check DeltaY
        if(is.null(iniList$DeltaY)){
            iniList$DeltaY <- iniList$oData$Header$DeltaXY[[2]]
            if(DL2)
                cat("'DeltaY' entry is set to", iniList$DeltaY, "\n")
        } else if (iniList$oData$Header$DeltaXY[[2]] != iniList$DeltaY){
            if(DL2)
                cat("Note: 'DeltaY' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        iniList$oData <- iniList$oData$Signal
        if(DL2) cat("\n")
    } else if (!(is.null(iniList$oData))){
        # check XSamples
        if (is.null(iniList$XSamples))
            iniList$XSamples <- nrow(iniList$oData)
        # check YSamples
        if (is.null(iniList$YSamples))
            iniList$YSamples <- ncol(iniList$oData)
        # check Xmin
        if (is.null(iniList$Xmin))
            iniList$Xmin <- -sqrt(0.5)*(ncol(iniList$rData)/iniList$XSamples)*0.5*(iniList$XSamples-1)
        # check Ymin
        if (is.null(iniList$Ymin))
            iniList$Ymin <- -sqrt(0.5)*(ncol(iniList$rData)/iniList$YSamples)*0.5*(iniList$YSamples-1)
        # check DeltaX
        if (is.null(iniList$DeltaX))
            iniList$DeltaX <- sqrt(0.5)*(ncol(iniList$rData)/iniList$XSamples)
        # check DeltaY
        if (is.null(iniList$DeltaY))
            iniList$DeltaY <- sqrt(0.5)*(ncol(iniList$rData)/iniList$YSamples)
    } else {
        if(is.null(iniList$XSamples))
            stop("'XSamples' entry missing in INI file")
        if(is.null(iniList$YSamples))
             stop("'YSamples' entry missing in INI file")
        if(is.null(iniList$Xmin))
             iniList$Xmin <- -sqrt(0.5)*(ncol(iniList$rData)/iniList$XSamples)*0.5*(iniList$XSamples-1)
        if(is.null(iniList$Ymin))
             iniList$Ymin <- -sqrt(0.5)*(ncol(iniList$rData)/iniList$YSamples)*0.5*(iniList$YSamples-1)
        if(is.null(iniList$DeltaX))
             iniList$DeltaX <- sqrt(0.5)*(ncol(iniList$rData)/iniList$XSamples)
        if(is.null(iniList$DeltaY))
             iniList$DeltaY <- sqrt(0.5)*(ncol(iniList$rData)/iniList$YSamples)
    }
    
} else {
    ###########################################################
    #  parameter checking for iterative reconstruction-method
    
    # load rData
    if (is.null(iniList$rData))
         stop("'rData' entry missing in INI file.")
    tmp <- unlist(strsplit(iniList$rData, ".", fixed=TRUE))
    fformatRadon <- tmp[length(tmp)]
    rfile <- iniList$rData
    if (DL2) cat("\n")
    iniList$rData <- readData(iniList$rData, DebugLevel)
    if (DL2) cat("\n")
    
    # If StartImage specified than loading StartImage
    if (iniList$StartImage != "None" ){
        tmp <- unlist(strsplit(iniList$StartImage, ".", fixed=TRUE))
        fformatStart <- tmp[length(tmp)]
        sfile <- iniList$StartImage
        iniList$StartImage <- readData(iniList$StartImage, DebugLevel)
        if (DL2) cat("\n")
    }

    ##############################################
    # checking the following parameters:
    # XSamples, YSamples, XMin, YMin, DeltaX and DeltaY
    if (iniList$StartImage != "None" && fformatStart %in% c("pet", "fif", "dat")){
        if(DL1) cat("Setting unknown parameters from '", sfile, "' \n", sep="")
        # check XSamples
        if(is.null(iniList$XSamples)){
            iniList$XSamples <- iniList$StartImage$Header$SignalDim[[1]]
            if(DL2)
                cat("'XSamples' entry is set to", iniList$XSamples, "\n")
        } else if (iniList$StartImage$Header$SignalDim[[1]] != iniList$XSamples){
            if(DL2)
                cat("Note: 'XSamples' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check YSamples
        if(is.null(iniList$YSamples)){
            iniList$YSamples <- iniList$StartImage$Header$SignalDim[[2]]
            if(DL2)
                cat("'YSamples' entry is set to", iniList$YSamples, "\n")
        } else if (iniList$StartImage$Header$SignalDim[[2]] != iniList$YSamples){
            if(DL2)
                cat("Note: 'YSamples' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check Xmin
        if(is.null(iniList$Xmin)){
            iniList$Xmin <- iniList$StartImage$Header$XYmin[[1]]
            if(DL2)
                cat("'Xmin' entry is set to", iniList$Xmin, "\n")
        } else if (iniList$StartImage$Header$XYmin[[1]] != iniList$Xmin){
            if(DL2)
                cat("Note: 'Xmin' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check Ymin
        if(is.null(iniList$Ymin)){
            iniList$Ymin <- iniList$StartImage$Header$XYmin[[2]]
            if(DL2)
                cat("'Ymin' entry is set to", iniList$Ymin, "\n")
        } else if (iniList$StartImage$Header$XYmin[[2]] != iniList$Ymin){
            if(DL2)
                cat("Note: 'Ymin' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check DeltaX
        if(is.null(iniList$DeltaX)){
            iniList$DeltaX <- iniList$StartImage$Header$DeltaXY[[1]]
            if(DL2)
                cat("'DeltaX' entry is set to", iniList$DeltaX, "\n")
        } else if (iniList$StartImage$Header$DeltaXY[[1]] != iniList$DeltaX){
            if(DL2)
                cat("Note: 'DeltaX' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        # check DeltaY
        if(is.null(iniList$DeltaY)){
            iniList$DeltaY <- iniList$StartImage$Header$DeltaXY[[2]]
            if(DL2)
                cat("'DeltaY' entry is set to", iniList$DeltaY, "\n")
        } else if (iniList$StartImage$Header$DeltaXY[[2]] != iniList$DeltaY){
            if(DL2)
                cat("Note: 'DeltaY' is different between '", fileName, "' and '", sfile, "'. INI-file is preferred. \n", sep="")
        }
        iniList$StartImage <- iniList$StartImage$Signal
        if(DL2) cat("\n")
    } else if (iniList$StartImage != "None"){
        # check XSamples
        if (is.null(iniList$XSamples))
            iniList$XSamples <- nrow(iniList$StartImage)
        # check YSamples
        if (is.null(iniList$YSamples))
            iniList$YSamples <- ncol(iniList$StartImage)
        # check Xmin
        if (is.null(iniList$Xmin))
            iniList$Xmin <- -0.5*((iniList$XSamples)-1)
        # check Ymin
        if (is.null(iniList$Ymin))
            iniList$Ymin <- -0.5*((iniList$YSamples)-1)
        # check DeltaX
        if (is.null(iniList$DeltaX))
            iniList$DeltaX <- 1
        # check DeltaY
        if (is.null(iniList$DeltaY))
            iniList$DeltaY <- 1
    }
    
    if(is.null(iniList$XSamples))
        stop("'XSamples' entry missing in INI file")
    if(is.null(iniList$YSamples))
        stop("'YSamples' entry missing in INI file")   
    
    ##############################################
    # checking the following parameters:
    # ThetaSamples, RhoSamples, ThetaMin, RhoMin, DeltaTheta and DeltaRho
    if (fformatRadon %in% c("pet", "fif", "dat")){
        if(DL1) cat("Setting unknown parameters from '", rfile, "' \n", sep="")
        # check ThetaSamples
        if(is.null(iniList$ThetaSamples)){
            iniList$ThetaSamples <- iniList$rData$Header$SignalDim[[1]]
            if(DL2)
                cat("'ThetaSamples' entry is set to", iniList$ThetaSamples, "\n")
        } else if (iniList$rData$Header$SignalDim[[1]] != iniList$ThetaSamples){
            if(DL2)
                cat("Note: 'ThetaSamples' is different between '", fileName, "' and '", rfile, "'. INI-file is preferred. \n", sep="")
        }
        # check RhoSamples
        if(is.null(iniList$RhoSamples)){
            iniList$RhoSamples <- iniList$rData$Header$SignalDim[[2]]
            if(DL2)
                cat("'RhoSamples' entry is set to", iniList$RhoSamples, "\n")
        } else if (iniList$rData$Header$SignalDim[[2]] != iniList$RhoSamples){
            if(DL2)
                cat("Note: 'RhoSamples' is different between '", fileName, "' and '", rfile, "'. INI-file is preferred. \n", sep="")
        }
        # check ThetaMin
        if(is.null(iniList$ThetaMin)){
            iniList$ThetaMin <- iniList$rData$Header$XYmin[[1]]
            if(DL2)
                cat("'ThetaMin' entry is set to", iniList$ThetaMin, "\n")
        } else if (iniList$rData$Header$XYmin[[1]] != iniList$ThetaMin){
            if(DL2)
                cat("Note: 'ThetaMin' is different between '", fileName, "' and '", rfile, "'. INI-file is preferred. \n", sep="")
        }
        # check RhoMin
        if(is.null(iniList$RhoMin)){
            iniList$RhoMin <- iniList$rData$Header$XYmin[[2]]
            if(DL2)
                cat("'RhoMin' entry is set to", iniList$RhoMin, "\n")
        } else if (iniList$rData$Header$XYmin[[2]] != iniList$RhoMin){
            if(DL2)
                cat("Note: 'RhoMin' is different between '", fileName, "' and '", rfile, "'. INI-file is preferred. \n", sep="")
        }
        # check DeltaTheta
        if(is.null(iniList$DeltaTheta)){
            iniList$DeltaTheta <- iniList$rData$Header$DeltaXY[[1]]
            if(DL2)
                cat("'DeltaTheta' entry is set to", iniList$DeltaTheta, "\n")
        } else if (iniList$rData$Header$DeltaXY[[1]] != iniList$DeltaTheta){
            if(DL2)
                cat("Note: 'DeltaTheta' is different between '", fileName, "' and '", rfile, "'. INI-file is preferred. \n", sep="")
        }
        # check DeltaRho
        if(is.null(iniList$DeltaRho)){
            iniList$DeltaRho <- iniList$rData$Header$DeltaXY[[2]]
            if(DL2)
                cat("'DeltaRho' entry is set to", iniList$DeltaRho, "\n")
        } else if (iniList$rData$Header$DeltaXY[[2]] != iniList$DeltaRho){
            if(DL2)
                cat("Note: 'DeltaRho' is different between '", fileName, "' and '", rfile, "'. INI-file is preferred. \n", sep="")
        }
        
        iniList$rData <- iniList$rData$Signal
        if(DL2) cat("\n")

    } else {
        # check ThetaSamples
        if (is.null(iniList$ThetaSamples))
            iniList$ThetaSamples <- nrow(iniList$rData)
        # check RhoSamples
        if (is.null(iniList$RhoSamples))
            iniList$RhoSamples <- ncol(iniList$rData)
        # check ThetaMin
        if (is.null(iniList$ThetaMin))
            iniList$ThetaMin <- 0
        # check RhoMin
        if (is.null(iniList$RhoMin))
            iniList$RhoMin <- -0.5*((2*round(sqrt(iniList$XSamples^2+iniList$YSamples^2)/2)+1)-1)
        # check DeltaTheta
        if (is.null(iniList$DeltaTheta))
            iniList$DeltaTheta <- pi/(iniList$ThetaSamples)
        # check DeltaRho
        if (is.null(iniList$DeltaRho))
            iniList$DeltaRho <- (2*abs(iniList$RhoMin)+1)/iniList$RhoSamples
    }

}
    
    if(is.null(iniList$XSamples))
        stop("'XSamples' entry missing in INI file")
    if(is.null(iniList$YSamples))
        stop("'YSamples' entry missing in INI file")
    if(is.null(iniList$Xmin))
        iniList$Xmin <- -0.5*((iniList$XSamples)-1)
    if(is.null(iniList$Ymin))
        iniList$Ymin <- -0.5*((iniList$YSamples)-1)
    if(is.null(iniList$DeltaX))
        iniList$DeltaX <- 1
    if(is.null(iniList$DeltaY))
        iniList$DeltaY <- 1
    
} else
    stop("mode='", mode,"' don't be supported. See help to 'iradon' or 'iradonIT' to get more information about mode.", sep="")
    

return(iniList)

}
