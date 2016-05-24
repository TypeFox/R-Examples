iniPetFifList <- function(listType="pet", data=NULL, imType="normal")
{
# ========================================================================
#
# iniPetFifList <- (listType, data)
#
# Description:
#     This function initialize a default header-list for writeData.
#
# listType  (character) Can be choosen between "pet" and "fif".
#           Set listType="pet" to generate a pet-list, with the following
#           structure:
#             Description  (optional) Can containing a short description of image.
#                          Maximal numbers of character has to be 80.
#             SignalDim    Containing the numbers of rows and columns of the 
#                          matrix. c(nrow(),ncol())
#             XYmin        Leftmost coor and lowest coor in original image.
#                          c(Xmin, Ymin)
#             DeltaXY      Quantisation steps in original image in X and Y 
#                          direction. c(DeltaX, DeltaY)
#           Set listType="fif" to generate a pet-list, with the following
#           structure:
#             FIFIdType    Id used to restore FIF:17737:'\0''\0''E''I'
#             FileName     Name used for saving/restoring this image.
#             Description  See above.
#             Date         (optional) Date (YYYY-MM-DD)
#             SignalDim    See above.
#             ArrayType    Defines number format. 1 for Real and 2 for Complex.
#             XYmin        See above.   
#             DeltaXY      See above.
#             SignalMinMax The lowest and highest signalvalue in matrix.
#                          c(min(),max())
#
# data      With it, it is possible to consign the matrix and to generate the
#           parameter SignalDim, XYmin, DeltaXY and SignalMinMax with default
#           values. In the other case, you have to specifiy these values.
#
#
#


  matrixTRUE <- 0
  if (!is.character(listType))
      stop("'listType' has to be of type character.")
  if ( !(is.null(data)) & !(is.matrix(data)) )
      warning("'data' has to be of type matrix. Default 'none' is used.")
  else if (is.matrix(data))
      matrixTRUE <- 1
  if (is.character(imType) & imType=="normal")
      imType <- TRUE
  else if(is.character(imType) & imType=="radon")
      imType <- FALSE
  else {
      warning("'imType' has to be 'normal' or 'radon'. Default 'normal' is used.")
      imType <- TRUE
  }

  if (listType == "pet"){
      iniList <- list()
      length(iniList) <- 4
      names(iniList) <- c("Description", "SignalDim", "XYmin", "DeltaXY")
      iniList$Description <- "" # OLD: "\0"
      if (matrixTRUE){
          iniList$SignalDim <- dim(data)
          if(imType){ 
              iniList$XYmin <- -0.5*(dim(data)-1)
              iniList$DeltaXY <- c(1,1)
          } else {
              iniList$XYmin <- c(0, -0.5*(ncol(data)-1))
              iniList$DeltaXY <- c(pi/(nrow(data)), 1)
          }
      }
        
  } else if (listType == "fif"){
      iniList <- list()
      length(iniList) <- 9
      names(iniList) <- c("FIFIdType", "FileName", "Description", "Date", 
          "SignalDim", "ArrayType", "XYmin", "DeltaXY", "SignalMinMax" )
    
      iniList$FIFIdType <- 17737
      iniList$FileName <- ""  # OLD: "\0"
      iniList$Description <- ""  # OLD: "\0"
      iniList$Date <- as.character(Sys.Date()) #format(Sys.Date(), "%d%m%y")
      iniList$ArrayType <- 1
      if (matrixTRUE){
          iniList$SignalDim <- dim(data)
          if(imType){ 
              iniList$XYmin <- -0.5*(dim(data)-1)
              iniList$DeltaXY <- c(1, 1)
          } else {
              iniList$XYmin <- c(0, -0.5*(ncol(data)-1))
              iniList$DeltaXY <- c(pi/(nrow(data)), 1)
          }
          iniList$SignalMinMax <- c(min(data), max(data))
      }
    
  } else
      stop("'listType' =",listType,"doesn't supported. \n")

  return(iniList)
}