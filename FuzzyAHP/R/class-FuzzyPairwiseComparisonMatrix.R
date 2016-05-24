#' Class "FuzzyPairwiseComparisonMatrix"
#' @description  An S4 class to represent a fuzzy pairwise comparison matrix.
#'
#' @slot fnMin A matrix of minimal values of fuzzy preferences.
#' @slot fnModal A matrix of modal values of fuzzy preferences.
#' @slot fnMax A matrix of maximal values of fuzzy preferences.
#' @slot variableNames Names of variables in the pariwise comparison matrix obtained either as colnames or rownames.
#'
#' @export
#' @include class-PairwiseComparisonMatrix.R
#'
setClass(
  Class="FuzzyPairwiseComparisonMatrix",

  representation(
    fnMin = "matrix",
    fnModal = "matrix",
    fnMax = "matrix",
    variableNames = "character"
  ),
  validity=function(object)
  {
    if(length(which(object@fnMin>object@fnModal)) != 0){
      return("Cannot create fuzzy pairwise comparison matrix. Minimal and modal values are not aligned correctly!")
    }

    if(length(which(object@fnModal>object@fnMax)) != 0){
      return("Cannot create fuzzy pairwise comparison matrix. Modal and maximal values are not aligned correctly!")
    }


  }
)





#'  Function that creates Fuzzy Pairwise Comparions Matrix
#'
#' @description
#' This methods construct object \linkS4class{FuzzyPairwiseComparisonMatrix} based on provided
#' \linkS4class{PairwiseComparisonMatrix} and an optional fuzzy scale.
#'
#' @details
#' Parameter fuzzyScale is expected as a vector containing n*3 values that represent triangular fuzzy
#' numbers used as fuzzy intensity of importance (only the values equal or higher than 1, inverse values
#' are calculated automatically). The values need to be ordered by fuzzy values. Default value of this
#' parameter is \code{as.double(c(1/2,1,2,1,2,3,2,3,4,3,4,5,4,5,6,5,6,7,6,7,8,7,8,9,8,9,9))}. Another
#' possibility is eg. \code{as.double(c(1/3,1,3,1,3,5,3,5,7,5,7,9,7,9,9))}
#'
#' If param \code{pairwiseComparisonMatrix} is \code{matrix} then it needs to be of \code{character} type.
#' Each element in the matrix must be specified as triplet \code{"x;y;z"}, where \code{x<=y<=z}.
#' From this matrix a pairwise comparison is constructed from \code{y} values and \code{x} and \code{z}
#' function as lower and upper limits of \code{y} respectively. In this case the optional parameter
#' \code{fuzzyScale} is not taken into account at all.
#'
#' @param pairwiseComparisonMatrix \linkS4class{PairwiseComparisonMatrix} or \code{matrix}.
#' @param fuzzyScale A numeric vector that definies fuzzy scale. Default scale is described in details.
#'
#' @usage fuzzyPairwiseComparisonMatrix(pairwiseComparisonMatrix, fuzzyScale)
#'
#' @return Object of class \linkS4class{FuzzyPairwiseComparisonMatrix}
#'
#' @export
#' @rdname fuzzyPairwiseComparisonMatrix-methods
#' @name fuzzyPairwiseComparisonMatrix
setGeneric("fuzzyPairwiseComparisonMatrix",
           signature = c("pairwiseComparisonMatrix"),
           function(pairwiseComparisonMatrix, fuzzyScale = as.double(c(1/2,1,2,
                                                                         1,2,3,
                                                                         2,3,4,
                                                                         3,4,5,
                                                                         4,5,6,
                                                                         5,6,7,
                                                                         6,7,8,
                                                                         7,8,9,
                                                                         8,9,9))
                    ) standardGeneric("fuzzyPairwiseComparisonMatrix"))

#' @rdname fuzzyPairwiseComparisonMatrix-methods
#' @aliases fuzzyPairwiseComparisonMatrix,PairwiseComparisonMatrix,fuzzyScale-method
setMethod(
  f="fuzzyPairwiseComparisonMatrix",
  signature(pairwiseComparisonMatrix = "PairwiseComparisonMatrix"),
  definition=function(pairwiseComparisonMatrix, fuzzyScale)
  {

    if(!(length(fuzzyScale)%%3==0)){
      stop(paste("The fuzzy scale lenght has to be x*3. This fuzzy scale does not fulfill this condition."))
    }

    size = nrow(pairwiseComparisonMatrix@values)
    # prepare 3 matrices for fuzzy values, each hase size of the original matrix
    fnMin = matrix(data = 0, nrow = size, ncol = size)
    fnModal = matrix(data = 0, nrow = size, ncol = size)
    fnMax = matrix(data = 0, nrow = size, ncol = size)

    # diagonal fuzzy values are crips ones
    v1diagonal = c(1,1,1)

    # transfor fuzzy scale into matrix
    #valuesMatrix = matrix(data = fuzzyScale, nrow = 9, ncol = 3, byrow = TRUE)
    valuesMatrix = matrix(data = fuzzyScale, nrow = length(fuzzyScale)/3, ncol = 3, byrow = TRUE)


    # prepare matrix of inverted fuzzy values from the scale
    invertedValuesMatrix = matrix(data = 1.0, nrow = length(fuzzyScale)/3, ncol = 3, byrow = TRUE ) / valuesMatrix
    invertedValuesMatrix = cbind(invertedValuesMatrix[,3], invertedValuesMatrix[,2], invertedValuesMatrix[,1])

    for (i in 1:size){
      for (j in 1:size){

        charIntensity = pairwiseComparisonMatrix@valuesChar[i,j]

        if (i==j & charIntensity == "1"){
          fnMin[i,j] = v1diagonal[1]
          fnModal[i,j] = v1diagonal[2]
          fnMax[i,j] = v1diagonal[3]
        }
        else if (nchar(charIntensity)==3 & substr(charIntensity,1,2)=="1/"){
          number = which(valuesMatrix[,2] == as.integer(substr(charIntensity,3,4))) #as.integer(substr(charIntensity,3,4))

          if(length(number)==0){
            stop(paste("Value ",charIntensity," does not exist in fuzzy scale!", sep = ""))
          }

          fnMin[i,j] = invertedValuesMatrix[number,1]
          fnModal[i,j] = invertedValuesMatrix[number,2]
          fnMax[i,j] = invertedValuesMatrix[number,3]
        }
        else if (nchar(charIntensity)==1){
          number = which(valuesMatrix[,2] == as.integer(substr(charIntensity,1,2))) #as.integer(substr(charIntensity,1,2))

          if(length(number)==0){
            stop(paste("Value ",charIntensity," does not exist in fuzzy scale!", sep = ""))
          }

          fnMin[i,j] = valuesMatrix[number,1]
          fnModal[i,j] = valuesMatrix[number,2]
          fnMax[i,j] = valuesMatrix[number,3]
        }
        else{
          stop("This should never happen. Error in format of pairwise comparison matrix.")
        }
      }
    }

    return(new("FuzzyPairwiseComparisonMatrix", fnMin = fnMin, fnModal = fnModal, fnMax = fnMax, variableNames = pairwiseComparisonMatrix@variableNames))
  }
)

#' @rdname fuzzyPairwiseComparisonMatrix-methods
#' @aliases fuzzyPairwiseComparisonMatrix,PairwiseComparisonMatrix,fuzzyScale-method
setMethod(
  f="fuzzyPairwiseComparisonMatrix",
  signature(pairwiseComparisonMatrix = "matrix"),
  definition=function(pairwiseComparisonMatrix, fuzzyScale)
  {
    if(typeof(pairwiseComparisonMatrix)!="character"){
      stop("Can only parse character matrix as fuzzy pairwise comparison matrix!")
    }

    if(nrow(pairwiseComparisonMatrix)!=ncol(pairwiseComparisonMatrix)){
      stop(paste("The fuzzy pairwise comparison matrix is not a square matrix. Dimensions are - ncol = ",
                   ncol(pairwiseComparisonMatrix), ", nrow = ", nrow(pairwiseComparisonMatrix), ".", sep = ""))
    }

    sep = ";"

    size = nrow(pairwiseComparisonMatrix)
    mMin = matrix(data = 0, nrow = size, ncol = size)
    mModal = matrix(data = 0, nrow = size, ncol = size)
    mMax = matrix(data = 0, nrow = size, ncol = size)

    rownames(mModal) = rownames(pairwiseComparisonMatrix)
    colnames(mModal) = colnames(pairwiseComparisonMatrix)

    for(i in 1:size){
      for(j in 1:size){

        string = pairwiseComparisonMatrix[i,j]
        stringSplit = strsplit(string, sep)

        if(length(stringSplit[[1]])!=3){
          stop(paste0("Element on position [", i, ",", j, "] ", string, " cannot be split into three values by separtor ",
                      sep, " !"))
        }else{
          vMin = eval(parse(text=stringSplit[[1]][1]))
          vModal = eval(parse(text=stringSplit[[1]][2]))
          vMax = eval(parse(text=stringSplit[[1]][3]))

          if(vMin>vModal || vModal>vMax){
            stop(paste0("Element on [",i,",",j,"] is not ordered correctly. ", vMin, "<=",
                       vModal, "<=", vMax))
          }

          mMin[i,j] = vMin
          mModal[i,j] = vModal
          mMax[i,j] = vMax
        }
      }
    }

    mModal = pairwiseComparisonMatrix(mModal)

    return(fuzzyPairwiseComparisonMatrix1(mMin, mModal, mMax))

    # return(new("FuzzyPairwiseComparisonMatrix", fnMin = mMin, fnModal = mModal@values,
    #            fnMax = mMax, variableNames = .getVariableNames(pairwiseComparisonMatrix)))
  }
)





#'  Function that creates Fuzzy Pairwise Comparions Matrix
#'
#' @description
#' This methods construct object \linkS4class{FuzzyPairwiseComparisonMatrix} based on provided
#' \linkS4class{PairwiseComparisonMatrix} and two matrices that form lower an upper significant values of the
#' \linkS4class{PairwiseComparisonMatrix} that form middle significant value.
#'
#' @details
#' This function allows user to specify fuzzy pairwise comparison matrix that is not based on fuzzy scale but
#' rely more on user's specification. The middle significant values have to be definied by
#' \linkS4class{PairwiseComparisonMatrix} to ensure some elementary properties. The significant values provided
#' to this function have to be correctly ordered and fuzzy numbers have to be reciprocal otherwise the function
#' fails.
#'
#' @param lowerValues A matrix of \code{"double"} that consists of lower significant values.
#' @param pairwiseComparisonMatrix \linkS4class{PairwiseComparisonMatrix} that consists of middle significant values.
#' @param upperValues A matrix of \code{"double"} that consists of upper significant values.
#'
#' @usage fuzzyPairwiseComparisonMatrix1(lowerValues, pairwiseComparisonMatrix, upperValues)
#'
#' @return Object of class \linkS4class{FuzzyPairwiseComparisonMatrix}
#'
#' @export
#' @rdname fuzzyPairwiseComparisonMatrix1-methods
#' @name fuzzyPairwiseComparisonMatrix1
setGeneric("fuzzyPairwiseComparisonMatrix1",
           function(lowerValues, pairwiseComparisonMatrix, upperValues)
           standardGeneric("fuzzyPairwiseComparisonMatrix1"))

#' @rdname fuzzyPairwiseComparisonMatrix1-methods
#' @aliases fuzzyPairwiseComparisonMatrix1,matrix,PairwiseComparisonMatrix,matrix-method
setMethod(
  f="fuzzyPairwiseComparisonMatrix1",
  signature(lowerValues = "matrix", pairwiseComparisonMatrix = "PairwiseComparisonMatrix", upperValues = "matrix"),
  definition=function(lowerValues, pairwiseComparisonMatrix, upperValues)
  {
    if (ncol(lowerValues)!=ncol(upperValues) || ncol(lowerValues)!=ncol(pairwiseComparisonMatrix@values)){
      stop("Matrices do not have the same number of columns!")
    }

    if (nrow(lowerValues)!=nrow(upperValues) || nrow(lowerValues)!=nrow(pairwiseComparisonMatrix@values)){
      stop("Matrices do not have the same number of rows!")
    }

    if(typeof(lowerValues)!="double" || typeof(upperValues)!="double"){
      stop("Matrices lower and upperValues need to be double type!")
    }

    size = nrow(pairwiseComparisonMatrix@values)

    for (i in 1:size){
      for (j in 1:size){

        if(!(lowerValues[i,j]<=pairwiseComparisonMatrix@values[i,j] && pairwiseComparisonMatrix@values[i,j]<=upperValues[i,j])){
          stop(paste("Element on [",i,",",j,"] is not ordered correctly. ",lowerValues[i,j], "<=",
                     pairwiseComparisonMatrix@values[i,j], "<=", upperValues[i,j], sep=""))
        }

        if(i == j && (lowerValues[i,j]!=1 || upperValues[i,j]!=1)){
          stop(paste("Elements on the main diagonal need to be (1,1,1). Element [",i,",",j,"] is (",
                     lowerValues[i,j],",",pairwiseComparisonMatrix@values[i,j],",", upperValues[i,j],").", sep = ""))
        }else if(lowerValues[i,j]!=(1/upperValues[j,i])){
          stop(paste("Elements [",i,",",j,"]  and [",j,",",i,"] are not reciprocal. ",lowerValues[i,j],
                     " is not reciprocal to ", upperValues[j,i],".", sep = ""))
        }
      }
    }

    return(new("FuzzyPairwiseComparisonMatrix", fnMin = lowerValues, fnModal = pairwiseComparisonMatrix@values,
               fnMax = upperValues, variableNames = pairwiseComparisonMatrix@variableNames))
  }
)
