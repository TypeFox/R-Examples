setClass (
  Class = "DiscernibilityMatrix",
  representation = representation(
    discernibilityMatrix = "array"
  ),
  validity = function(object){
    #cat("~~~ DiscernibilityMatrix: inspector ~~~ \n")
    if(length(dim(object@discernibilityMatrix)) != 3){
      stop ("[DiscernibilityMatrix: validation] Discernibility Matrix must be an 3 dim array")      
    }else{}
    return(TRUE)
  }
)
#*******************************************************
#CONSTRUCTOR
setMethod (
  f="initialize",
  signature="DiscernibilityMatrix",
  definition=function(.Object,discernibilityMatrix){
    #cat ("~~~~~ DiscernibilityMatrix: initializator ~~~~~ \n")
    
    if(!missing(discernibilityMatrix)){
      #dimnames are added by DecisionTable
      .Object@discernibilityMatrix <- discernibilityMatrix
      validObject(.Object)# call of the inspector
    }else{
      .Object@discernibilityMatrix <- matrix(nrow=0,ncol=0)
    }
    return(.Object)
  }
)

#CONSTRUCTOR (USER FRIENDLY)
discernibilityMatrix <- function(theDiscernibilityMatrix){
  #cat ("~~~~~ DiscernibilityMatrix: constructor ~~~~~ \n")
  new (Class="DiscernibilityMatrix",discernibilityMatrix = theDiscernibilityMatrix)
}

#*******************************************************
#ACCESSORS

#Returns the complete decision table as array
setGeneric("getDiscernibilityMatrix",function(object){standardGeneric ("getDiscernibilityMatrix")})
setMethod("getDiscernibilityMatrix","DiscernibilityMatrix",
          function(object){
            return(object@discernibilityMatrix)
          }
)

#*******************************************************
#GENERIC METODS

#Summary

setMethod ("print","DiscernibilityMatrix",
           function(x,...){
             cat("*** Class DiscernibilityMatrix, method Print *** \n")
             disMat <- x@discernibilityMatrix
             if(length(disMat) != 0){
               mat <- .buildPrintMatrix(disMat)
               print(formatC(mat),quote=FALSE)
             }else{}
             cat("******* End Print (DiscernibilityMatrix) ******* \n")
           }
)


setMethod("show","DiscernibilityMatrix",
          function(object){
            cat("*** Class DiscernibilityMatrix, method Show *** \n")
            cat("* DiscernibilityMatrix (limited to a matrix 10x10) = \n")
            disMat <- object@discernibilityMatrix
            if(length(disMat) != 0){
              mat <- .buildPrintMatrix(disMat)
              nrowShow <- min(10,nrow(mat))
              ncolShow <- min(10,ncol(mat))
              print(formatC(mat[1:nrowShow,1:ncolShow]),quote=FALSE)
            }else{}
            cat("******* End Show (DiscernibilityMatrix) ******* \n")
          }
)


#*******************************************************
#METHODS

### computeCore
### - theArray is a 3 dim array representing the discernibility matrix of a decision table
###
###
### computeCore  returns a numeric vector indicating the conditions ids (columns ids) which are the core of the Decision Table object from which the Discernibility Matrix object was created
###
setGeneric (name = "computeCore",def = function(object){standardGeneric("computeCore")})
setMethod(
  f = "computeCore",
  signature = "DiscernibilityMatrix",
  definition = function(object){
    dm <- object@discernibilityMatrix
    core <- .computeCore(dm)
    return(core)
  }
)


#*******************************************************
#UTIL

### .buildPrintMatrix
### - theArray is a 3 dim array representing the discernibility matrix of a decision table
###
###
### .buildPrintMatrix returns a formated caracter matrix for printing and showing
###
.buildPrintMatrix <- function(theArray){
  disMat <- theArray
  lenI <- dim(disMat)[1]
  lenJ <- dim(disMat)[2]
  lenK <- dim(disMat)[3]
  if(length(disMat) != 0){
    mat <- matrix(NA,nrow = lenI, ncol = lenJ)
    for(i in 1:lenI){
      for(j in 1:lenJ){
        ij <- disMat[i,j,]
        partialElement <- vector(mode = "numeric",length = 0)
        partialCounter <- 0
        for(k in 1:lenK){
          if(!is.na(ij[k])){
            if(ij[k]){
              partialCounter <- (partialCounter +1)
              partialElement[partialCounter] <- k
            }else{}
          }else{}
        }
        if(length(partialElement) > 0){
          mat[i,j] <- paste("C",partialElement,sep="",collapse=",")
        }else{}
      }
    }
    mat <- t(mat)#Traspones juts to make it look nicer
    rownames(mat) <- paste("R",1:nrow(mat),sep="")
    colnames(mat) <- paste("R",1:nrow(mat),sep="")
  }else{}  
  return(mat)
}


### .computeCore
### - theDiscernibilityMatrix is a 3D array representing all rule's condition differences of a decision table
###
###
### .computeCore returns a numeric vector  indicating the decision table columns id which are part of the core
###
.computeCore <- function(theDiscernibilityMatrix){
  ruleCount <- nrow(theDiscernibilityMatrix)
  discernibilityMatrix <- theDiscernibilityMatrix
  
  #********************************************************************************************************************************************  
  #PARALELL COMPUTING - SNOWFALL
  #********************************************************************************************************************************************
  idRedMat <- apply(discernibilityMatrix,c(1,2),.reductIdentification)#The result is a matrix, it neeeds processing
  #********************************************************************************************************************************************
  
  idRedVec <- vector(mode = "numeric", length = 0)
  
  for(i in 1:ruleCount){
    for(j in 1:ruleCount){
      if(!is.na(idRedMat[i,j])){
        idRedVec[length(idRedVec) + 1] <- idRedMat[i,j]
      }else{}
    }
  }
  idRedVec <- unique(idRedVec)
  idRedVec <- sort(idRedVec)
  return(idRedVec)
}


### .reductIdentification
### - discernibilityMatrixElement is a boolean vector representing an element of the Discernibility Matrix
###
###
### .reductIdentification returns a boolean vector indicating the position of the condition which is part of the core
###
.reductIdentification <- function(discernibilityMatrixElement){
  conditionCount <- length(discernibilityMatrixElement)
  result <- NA
  #The condition being part of the reduct is identified cause it is alone in the discernibilityMatrixElement
  if(sum(discernibilityMatrixElement,na.rm = TRUE) == 1){
    #Identifies the condition which is part of the reduct
    for(i in 1:conditionCount){
      if(discernibilityMatrixElement[i] == TRUE){
        result <- i
        break
      }else{}
    }
  }else if(sum(discernibilityMatrixElement,na.rm = TRUE) != 1){
    result <- NA
  }else{}
  return(result)
}