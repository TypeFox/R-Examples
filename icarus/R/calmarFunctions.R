# copyright (C) 2014-2016 A.Rebecq
# Functions designed so that calibration can be made in a familiar
# setting for Calmar and Calmar2 users

# Remarque : calmarMatrix = "matrice des marges sans la colonne des noms"

nModalities = function(col)
{
  return(length(unique(col)))
}

calibrationMatrix = function(entryMatrix, popVector=TRUE, isQuantitative=NULL)
{
  if(is.null(isQuantitative)) {
    isQuantitative <- rep(FALSE, ncol(entryMatrix))
  }

  entryMatrix = data.matrix(entryMatrix)

  # Initialization of return matrix
  nRows = nrow(entryMatrix)
  nCols = 0

  N = ncol(entryMatrix)

  # Particular case if entryMatrix has only one row
  if(is.null(N)) {
    N=1
    nRows = length(entryMatrix)
  }

  for(i in 1:N)
  {
    nCols = nCols + nModalities(entryMatrix[,i])
  }

  namesMatrix = names(entryMatrix)
  calibrationMatrix = matrix(0, nRows, 0, byrow=T)
  for(i in 1:N)
  {
    if(!isQuantitative[i]) {
      calibrationMatrix = cbind(calibrationMatrix, colToDummies(entryMatrix[,i], namesMatrix[i]))
    } else {
      calibrationMatrix = cbind(calibrationMatrix, entryMatrix[,i])
    }
  }

  # Add "population" vector
  if(popVector) {
    calibrationMatrix = cbind(calibrationMatrix, rep(1,nrow(calibrationMatrix)))
  }

  return(calibrationMatrix)
}

dummyModalitiesMatrix = function(entryMatrix)
{
  dmatrix = calibrationMatrix(entryMatrix)
  dmatrix[dmatrix!=0] = 1
  return(dmatrix)
}

# TODO : move out of calmarFunctions ?
# (or at least should be "private")
HTtotals = function(dummyModalitiesMatrix, weights)
{
  return(weights%*%dummyModalitiesMatrix)
}

# "createCalibrationMatrix" ensures compatibility with first version of icarus 
# (then called gaston 0.0.1)
createCalibrationMatrix = function(marginMatrix, data, popVector=TRUE)
{
  # Selection des variables de calage dans la table
  # (ainsi que leur caractere qualitatif / quantitatif)
  selectVector = marginMatrix[,1]
  isQuantitative = as.numeric(marginMatrix[,2])

  isQuantitative[isQuantitative != 0] <- 1
  isQuantitative <- 1 - as.numeric(isQuantitative) # is considered as boolean by R

  Xs = data[,selectVector]
  # Mise en forme de la matrice de calage
  matrixCal = calibrationMatrix(Xs, popVector, isQuantitative)

  return(matrixCal)
}

formatMargins = function(calmarMatrix, calibrationMatrix, popTotal=NULL, pct=FALSE)
{
  # Create empty vector of margins
  cMatrixCopy = calmarMatrix
  if(is.vector(cMatrixCopy)) {
    cMatrixCopy = t(as.matrix(calmarMatrix))
    calmarMatrix = t(as.matrix(calmarMatrix))
  }
  typeMargins = cMatrixCopy[,1]
  typeMargins[typeMargins==0] = 1
  cMargins = rep(0,sum(typeMargins))

  # Fill cMargins
  i=1
  curRow = 1
  while(curRow <= nrow(calmarMatrix))
  {
    if(calmarMatrix[curRow,1] == 0)
    {
      cMargins[i]=calmarMatrix[curRow,2]
      i=i+1
    }
    else
    {
      ## TODO : change by using parameter pct
      ## ... this means that pct is considered TRUE by default
#       if(!is.null(popTotal)) {
#         popTotalNum <- popTotal
#       } else {
#         popTotalNum <- 1
#       }

      n = calmarMatrix[curRow,1]

      ## TODO : change to pct
      ## If categorial margins are not entered as percentages,
      ## do not multiply by popTotal (except if it is popVector !)
      
      if( all(calmarMatrix[curRow,2:(n+1)] < 1) && (is.null(popTotal) || !pct) ) {
        warning(paste("All margins in variable ",curRow,"are less than 1 : should they be considered as percentages ?"))
      }
      
#       if( all(calmarMatrix[curRow,2:(n+1)] >= 1) ) {
#         popTotalNum <- 1
#       }
      
      if(pct) {
        if(is.null(popTotal)) {
          stop("popTotal has to be set when pct is TRUE")
        } else {
          popTotalNum <- popTotal
        }
      } else {
        popTotalNum <- 1
      }
      
      for(j in 2:(n+1))
      {
        cMargins[i] = calmarMatrix[curRow,j]*popTotalNum
        i = i+1
      }
    }
    curRow = curRow+1
  }

  # If there is still one column, it is the population one, so we add popTotal to cMargins
  # ... unless specified otherwise
  if(i <= ncol(calibrationMatrix) && !is.null(popTotal))
    cMargins[i] = popTotal

  return(cMargins)
}

#' Stats for initial weights, calibrated weights, and margins.
#' @description
#' Gives stats about the calibration process: differences between 
#' totals after/before calibration and margins. Totals for categorical
#' variables are displayed in percentages.
#' (same as first panels output in Calmar/Calmar 2)
#' Output is a list, which might not be convenient for exports (e.g. for integration
#' into a scientific report). In such cases,
#' use function \code{\link{marginStats}}, which outputs a dataframe.
#' @param data dataframe containing the survey data
#' @param marginMatrix matrix of margins
#' @param popTotal total of population, useful if margins are entered in relative value
#' @param pct Set this to true if margins for categorical variables are written in percentages
#' @param colWeights name of weights column in the dataframe
#' @param colCalibratedWeights name of calibrated weights column in the dataframe (if applicable)
#' @param calibThreshold If difference between calibration estimate and margin differ more than
#' this parameter, calibration is considered to have failed
#' @return List containing stats on weights and margins
#' @seealso \code{\link{marginStats}}
#' @keywords statistics, stats, description, results
#' @export
calibrationMarginStats = function(data, marginMatrix, popTotal=NULL, pct=FALSE, colWeights, colCalibratedWeights=NULL, calibThreshold=1.0) {

  displayCalibratedWeights <- TRUE

  if(is.null(colCalibratedWeights)) {
    displayCalibratedWeights <- FALSE
    colCalibratedWeights <- colWeights
  }

  if(displayCalibratedWeights) {
    textAfter <- "After Calibration"
  } else {
    textAfter <- "Current"
  }

  enteredAsPct <- FALSE
  popTotalMarginDisplay <- popTotal
  if(is.null(popTotal)) {
    enteredAsPct <- FALSE
    if(displayCalibratedWeights) {
      popTotal <- sum(data[colCalibratedWeights])
    } else {
      popTotal <- sum(data[colWeights])
    }
    popTotalMarginDisplay <- NA
  }
  
  if(pct) {
    enteredAsPct <- TRUE
  }

  toWarn = FALSE
  displayWarningMessage = FALSE

  # Somme des poids (total)
  totalWeights = sum(data.matrix(data[colWeights]))
  totalCalibrated = sum(data[colCalibratedWeights])

  vecTotal = c(totalWeights, totalCalibrated, popTotalMarginDisplay)
  names(vecTotal) = c("Before calibration",textAfter, "Margin")

  vecTotal = round(vecTotal,2)

  marginStatsList = list(vecTotal)

  marginNames = marginMatrix[,1]

  if(is.null(marginMatrix)) {
    names(marginStatsList) = c("Total")
    return(marginStatsList)
  }

  # Other margins
  for(i in 1:nrow(marginMatrix)) {

    toWarn = FALSE
    vecTotal = NULL

    if(as.numeric(marginMatrix[i,2]) == 0) { # If variable is numeric

      sumWeights = data.matrix(data[marginNames[i]])[,1] %*% data.matrix(data[colWeights])[,1]
      sumCalibrated = data.matrix(data[marginNames[i]])[,1] %*% data.matrix(data[colCalibratedWeights])[,1]
      margin = as.numeric(marginMatrix[i,3])

      vecTotal = c(sumWeights, sumCalibrated, margin)
      vecTotal = as.numeric(vecTotal)
      vecTotal = round(vecTotal,2)

      # Check if calibration is exact
      if(is.na(sumCalibrated)) stop(paste("Modality is present in margin tables but not in sample : ",i,";",j))
      if(abs(sumCalibrated - margin) >= calibThreshold) {
        toWarn = TRUE
        displayWarningMessage = TRUE
        #vecTotal = c(vecTotal,"*") # Old convention (same as in Calmar)
        vecTotal = c(vecTotal,round(abs((sumCalibrated - margin)/margin),4))
      }

      if(toWarn == FALSE) {
        names(vecTotal) = c("Before calibration",textAfter,"Margin")
      } else {
        names(vecTotal) = c("Before calibration",textAfter,"Margin", "Warning")
      }
    } else { # If variable has modalities
      modalities = data.matrix(unique(data[marginNames[i]])[,1])
      modalities = sort(modalities)

      # TODO : Assert length(modalities) == marginMatrix[i,2]

      for(j in 1:marginMatrix[i,2]) {

        toWarn = FALSE
        sumWeights = sum(data.matrix(data[data[marginNames[i]] == modalities[j],][colWeights]))
        sumCalibrated = sum(data.matrix(data[data[marginNames[i]] == modalities[j],][colCalibratedWeights]))

        if(!enteredAsPct) {
          ## By convention, margin for categorical variables are given in percentages
          margin = as.numeric(marginMatrix[i,2+j])
          # tempStatVec = c(sumWeights, sumCalibrated, margin)
          tempStatVec = c(sumWeights/totalWeights*100, sumCalibrated/totalCalibrated*100, margin/popTotal*100)
        } else {
          margin = as.numeric(marginMatrix[i,2+j])
          tempStatVec = c(sumWeights/totalWeights*100, sumCalibrated/totalCalibrated*100, margin*100)
        }

        #tempStatVec = c(sumWeights, sumCalibrated, margin) # TODO : change here level / structure

        # tempStatVec = c(sumWeights/totalWeights*100, sumCalibrated/totalCalibrated*100, margin/popTotal*100)

        tempStatVec = round(tempStatVec,2)


        # Check if calibration is exact
        if(is.na(sumCalibrated)) stop(paste("Modality is present in margin tables but not in sample : ",i,";",j))
        if(abs(sumCalibrated - margin) >= calibThreshold) {
#           toWarn = TRUE
          displayWarningMessage = TRUE
#           tempStatVec = c(tempStatVec, "*")
        }

        vecTotal = rbind(vecTotal, tempStatVec, deparse.level = 0)
      }

      # rownames = marginName_modalities(i)
      rownames(vecTotal) = modalities

      # "Little stars" if not perfectly calibrated
      if(toWarn == FALSE) {
        colnames(vecTotal) = c("Before calibration",textAfter,"Margin")
      } else {
        colnames(vecTotal) = c("Before calibration",textAfter,"Margin", "Warning")
      }



    }

    marginStatsList[[i+1]] = vecTotal
  }


  # Name of statsMargesList
  names(marginStatsList) = c("Total", marginNames)

  if(displayWarningMessage && displayCalibratedWeights)
    writeLines("Careful, calibration may not be exact")

  return(marginStatsList)
}

#' Stats for initial weights, calibrated weights, and margins.
#' @description
#' Just like \code{\link{calibrationMarginStats}}, gives stats about the calibration process: 
#' differences between totals after/before calibration and margins. Totals for categorical
#' variables are displayed in percentages. The last column, named "difference", shows
#' the difference (in percentage points) between initial estimates and margins (if colCalibratedWeights is NULL) 
#' or between calibrated estimates and margins (if colCalibratedWeights is not NULL).
#' Output is a dataframe, which might be more convenient to export than a list
#' (e.g. for integration into reports).
#' @param data dataframe containing the survey data
#' @param marginMatrix matrix of margins
#' @param pct Set this to true if margins for categorical variables are written in percentages
#' @param popTotal total of population, useful if margins are entered in relative value
#' @param colWeights name of weights column in the dataframe
#' @param colCalibratedWeights name of calibrated weights column in the dataframe (if applicable)
#' @param calibThreshold If difference between calibration estimate and margin differ more than
#' this parameter, calibration is considered to have failed
#' @return Dataframe containing stats on weights and margins
#' @seealso \code{\link{calibrationMarginStats}}
#' @keywords statistics, stats, description, results
#' @export
marginStats <- function(data, marginMatrix, pct=FALSE, popTotal=NULL, colWeights
                        , colCalibratedWeights=NULL, calibThreshold=1.0) {

  listMarginStats <- calibrationMarginStats(data, marginMatrix, popTotal, pct, colWeights
                                            , colCalibratedWeights, calibThreshold)
  marginStatsDF <- do.call(rbind.data.frame, listMarginStats)

  ## Compute column difference
  marginStatsDF <- marginStatsDF[,-c(4)]
  if( is.null(colCalibratedWeights) ) {

    marginStatsDF <- marginStatsDF[,-c(2)] # Do not display calibrated weigths column
    
    marginStatsDF[,3] <- round(abs(marginStatsDF[,2] - marginStatsDF[,1])/marginStatsDF[,2]*100,2)
    
    ## Correct coefficients for categorical variables
    marginStatsDF <- correctCoefsCategorical(marginStatsDF, marginMatrix)
    
    names(marginStatsDF) <- c("Before calibration","Margin", "Difference (pct)")
    
    
  } else {
    
    marginStatsDF[,4] <- round(abs(marginStatsDF[,3] - marginStatsDF[,2])/marginStatsDF[,3]*100,2)
    
    ## Correct coefficients for categorical variables
    marginStatsDF <- correctCoefsCategorical(marginStatsDF, marginMatrix, ncol1=2, ncol2=3, ncol3=4)
    
    colnames(marginStatsDF) <- c("Before calibration","After calibration","Margin","Difference (pct)")
  }

  return(marginStatsDF)

}

## Private function, used in marginMatrix to account for
## categorical variables, whose stats are displayed in percentages
correctCoefsCategorical <- function(marginStatsDF_init, marginMatrix, ncol1=1, ncol2=2, ncol3=3) {
  
  marginStatsDF <- marginStatsDF_init
  
  nModalCateg <- 0
  for(i in 1:nrow(marginMatrix)) {
    nModal <- as.numeric(marginMatrix[i,2]) 
    if(nModal > 0) {
      
      for(j in 1:(nModal)) {
        ## Offset of 1 because of popTotal in first line of marginStatsDF
        marginStatsDF[i+nModalCateg+1,ncol3] <- round(abs(marginStatsDF[i+nModalCateg+1,ncol2] - marginStatsDF[i+nModalCateg+1,ncol1]),2)
        if(j < nModal) nModalCateg <- nModalCateg + 1
      }
    }
  }
  
  return(marginStatsDF)
  
}

## TODO : deprecate, never used
## Replaced by checkNumberMargins
# Check validity of marginMatrix
checkMarginMatrix = function(marginMatrix) {

  checkMatrix = FALSE

  if(is.null(marginMatrix)) return(TRUE) # Case NULL is OK

  # TODO :
  # Check if there are : 1 names column, 1 modalities column and
  # n other columns with n = max(modalities)


  # Check if sum(lines where modalities >=2) = 1.000000


  return(checkMatrix)
}

# Displays number of NAs among margins
missingValuesMargins = function(data, marginMatrix) {

  nVar = nrow(marginMatrix)
  marginNames = marginMatrix[,1]
  returnMatrix = cbind(marginNames, rep(0,nVar))

  for(i in 1:nVar) {
    returnMatrix[i,2] = nrow(data[is.na(data[marginNames[i]]),])
  }

  colnames(returnMatrix) = c("Margin","Missing values")

  return(returnMatrix)
}

# Checks if number of modalities in data matches expected ones according
# to marginMatrix
checkNumberMargins = function(data, marginMatrix) {

  returnBool = TRUE
  marginNames = marginMatrix[,1]

  for(i in 1:length(marginNames)) {

    nModalities = length(table(data.matrix(data[marginNames[i]])))
    expectedModalities = as.numeric(marginMatrix[i,2])
    if(nModalities != expectedModalities && expectedModalities > 0) { ## "0" indicates calibration is made on quantitative total
      writeLines(paste("Error on column ",marginNames[i]," : ",nModalities," modalities in data and ",expectedModalities," expected in margins"))
      return(FALSE)
    }

  }

  return(TRUE)
}

#' Regroup calibration modalities
#' @description 
#' Beware, this function modifies the calibrationMatrix and marginMatrix objects entered in parameter?
#' Regroups modalities entered in "vecModalities" into single
#' "newModality" in "calibrationMatrix" and adapts "marginMatrix" to the new concept.
#' Typical usage is right before a calibration (and after comptutation of marginMatrix), when
#' you realise calibration output is better when several modalities are reduced to one.
#' (typically very rare modalities, on which calibration constraints are very restrictive).
#' Uses pseudo-"call by reference" via eval.parent because 2 objects are modified :
#' calibrationMatrix and marginMatrix
#' @param calibrationMatrix calibration matrix
#' @param marginMatrix matrix containing the margins to the Icarus format
#' @param calibrationVariable name of the calibration varaible for which regroupment has to be done
#' @param vecModalities Initial modalities of the variable
#' @param newModality Regrouped modalities of the variable
#' 
#' @examples
#' ## Suppose we have a calibration matrix and a margin matrix containing information
#' ## for two categorical variables "X1" (10 modalities) and "X2" (5 modalities)
#' 
#' matrixCal <- data.frame(matrix(
#'                c(floor(10*runif(100))+1,floor((5)*runif(100))+1,
#'                floor(10*runif(100))+1,rep(10,100)),
#'                ncol=4))
#' marginMatrix <- matrix(c("X1",10,rep(1/10,10),
#'                  "X2",5,rep(1/5,5),rep(0,5)), nrow=2, byrow=TRUE)
#' 
#' # table(matrixCal$X1)
#' # 1  2  3  4  5  6  7  8  9 10 
#' # 9  8  8  8 11 15 13  6 10 12 
#' # marginMatrix
#' # [,1] [,2] [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]  [,10] [,11] [,12]
#' # [1,] "X1" "10" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1"
#' # [2,] "X2" "5"  "0.2" "0.2" "0.2" "0.2" "0.2" "0"   "0"   "0"   "0"   "0" 
#' 
#' regroupCalibrationModalities(matrixCal, marginMatrix, "X1", c(3,4,8), "0")
#' 
#' # table(matrixCal$X1)
#' # 0  1  2  5  6  7  9 10 
#' # 22  9  8 11 15 13 10 12 
#' # marginMatrix
#' # [,1] [,2] [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]  [,10]
#' # [1,] "X1" "8"  "0.3" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1" "0.1"
#' # [2,] "X2" "5"  "0.2" "0.2" "0.2" "0.2" "0.2" "0"   "0"   "0"  
#' 
#' @export
regroupCalibrationModalities <- function(calibrationMatrix, marginMatrix, calibrationVariable, vecModalities, newModality) {
  
  # First, check if number of modalities match in calibrationMatrix and marginMatrix,
  # otherwise stop
  if(!checkNumberMargins(calibrationMatrix, marginMatrix))
    stop("Number of modalities must match between calibrationMatrix and marginMatrix to regroup calibration modalities.")
  
  newCalibrationMatrix <- calibrationMatrix
  newMarginMatrix <- marginMatrix
  
  ## Modification in calibrationMatrix
  newCalibrationMatrix[calibrationVariable] <- regroupUnContiguuousModalities(data.matrix(newCalibrationMatrix[calibrationVariable]), vecModalities, newModality)
  
  ## Modification in marginMatrix
  calVarModalities <- unique(data.matrix(calibrationMatrix[calibrationVariable]))
  
  if(newModality %in% calVarModalities) {
    stop("New modality cannot be a modality that already exists in calibration matrix")
  }
  
  orderedCalVarModalities <- calVarModalities[order(calVarModalities)]
  indicesVecModalities <- which(orderedCalVarModalities %in% vecModalities)
  indicesVecModalities <- indicesVecModalities+2 ## First two columns are name and nModalities
  
  modifiedLine <- marginMatrix[marginMatrix[,1] == calibrationVariable,]
  sumRegrouped <- sum(as.numeric(modifiedLine[indicesVecModalities]))
  modifiedLine <- modifiedLine[-indicesVecModalities]
  
  # Insert new margin (sum) to the right place
  modifiedLine <- modifiedLine[modifiedLine != 0]
  newCalVarModalities <- unique(data.matrix(newCalibrationMatrix[calibrationVariable]))
  orderedNewCalVarModalities <- newCalVarModalities[order(newCalVarModalities)]
  
  insertPosition <- which(orderedNewCalVarModalities==newModality)
  modifiedLine <- c(modifiedLine[1:(2+insertPosition-1)],sumRegrouped,
                    modifiedLine[(2+insertPosition):length(modifiedLine)])
  
  newNModalities <- as.numeric(modifiedLine[2]) - length(vecModalities) + 1
  modifiedLine[2] <- newNModalities
  
  
  # Add 0s to end line
  modifiedLine <- modifiedLine[1:(as.numeric(modifiedLine[2])+2)]
  modifiedLine <- c(modifiedLine, rep("0.0000",ncol(marginMatrix) - length(modifiedLine)))
  
  # Careful, sum of weights must be equal to 1 even after modalities have been regrouped
  sumMarginLine <- sum(as.numeric(modifiedLine[3:length(modifiedLine)]))
  
  if( sumMarginLine != 1 ) {

    maxMarginValue <- max(as.numeric(modifiedLine[3:(as.numeric(modifiedLine[2])+2)]))
    maxIndex <- which(as.numeric(modifiedLine[3:length(modifiedLine)]) == maxMarginValue)
    modifiedLine[maxIndex+2] <- maxMarginValue + 1 - sumMarginLine
  }
  
  # Replace in marginMatrix
  newMarginMatrix[marginMatrix[,1] == calibrationVariable,] <- modifiedLine
  
  # Check if last column of margin matrix is all 0s. If it is, drop last column
  # (means larger line has been reduced). Continue to do so until last colmun is not only 0s.
  while( sum(as.numeric(newMarginMatrix[,ncol(newMarginMatrix)])) == 0 ) {
    newMarginMatrix <- newMarginMatrix[, -ncol(newMarginMatrix)]
  }
  
  eval.parent(substitute(calibrationMatrix <- newCalibrationMatrix))
  eval.parent(substitute(marginMatrix <- newMarginMatrix))
}


#' Adds a margin to marginMatrix
#'
#' @param marginMatrix The matrix of margins to add the new margin to
#' @param varName Name of variable in calibration matrix corresponding
#' to the new margin
#' @param vecTotals values of margins (Calmar style) for the variable.
#' Note : if length(vecTotals) > 1, then sum(thresholdAdjustToOne) has to be 1.
#' @param adjustToOne if TRUE and sum(vecTotals) is nearly 1, modify values of vecTotals
#' so that sum is 1.
#' @param thresholdAdjustToOne adjust sum(vecTotals) to 1 if difference
#' is under thresholdAdjustToOne
#'
#' @export
addMargin <- function(marginMatrix, varName, vecTotals, adjustToOne=TRUE, thresholdAdjustToOne = 0.01) {

  if(varName %in% marginMatrix[,1]) {
    stop(paste(varName,"is already in margin matrix."))
  }
  
  newMarginMatrix <- marginMatrix

  # Length of vecTotals :
  if( length(vecTotals) == 1 ) {
    nModality <- 0
  } else {
    if( length(vecTotals) > 1 ) {
      nModality <- length(vecTotals)
    } else {
      stop("vecTotals must be non NULL vector")
    }
  }

  # TODO : adjust vecTotals to 1
  if( nModality > 1 && sum(vecTotals) != 1  ) {

    if(adjustToOne && abs(sum(vecTotals) - 1) < thresholdAdjustToOne) {
      # TODO : adjust highest value
      maxMarginValue <- max(as.numeric(vecTotals))
      maxIndex <- which.max(as.numeric(vecTotals))
      vecTotals[maxIndex] <- maxMarginValue + 1 - sum(vecTotals)
    } else {
      stop("sum(vecTotals) must be equal to 1.")
    }

  }

  newMarginLine <- c(varName, nModality, vecTotals)

  # newMarginLine must have right format before it is added to
  # newMarginMatrix
  if(length(newMarginLine) < ncol(newMarginMatrix)) {
    # Add missing zeroes :
    missingZeroes <- rep(0, ncol(newMarginMatrix) - length(newMarginLine))
    newMarginLine <- c(newMarginLine, missingZeroes)
  }

  if(length(newMarginLine) > ncol(newMarginMatrix)) {
    # Add columns of 0s to newMarginMatrix
    missingZeroes <- matrix(nrow = nrow(newMarginMatrix), ncol = (length(newMarginLine) - ncol(newMarginMatrix)), 0)
    newMarginMatrix <- cbind(newMarginMatrix, missingZeroes)
  }

  # Append to newMarginMatrix :
  newMarginMatrix <- rbind(newMarginMatrix, newMarginLine, deparse.level = 0)

  return(newMarginMatrix)
}

## Modifies margin
modifyMargin <- function(marginMatrix, varName, vecTotals, adjustToOne=TRUE, thresholdAdjustToOne = 0.01) {

  # Delete selected margin
  indexSelectedMargin <- NULL
  i <- 1
  while(i <= nrow(marginMatrix)) {
    if(marginMatrix[i,1] == varName) {
      indexSelectedMargin <- i
    }
    i <- i+1
  }

  newMarginMatrix <- marginMatrix[-indexSelectedMargin,]
  if(is.null(ncol(newMarginMatrix))) {
    newMarginMatrix <- t(as.matrix(newMarginMatrix))
  }

  # Add selected margin
  newMarginMatrix <- addMargin(newMarginMatrix, varName, vecTotals, adjustToOne, thresholdAdjustToOne)

  return(newMarginMatrix)
}

## Private function that creates margins to the right format
createFormattedMargins <- function(data, marginMatrix, popTotal=NULL, pct=FALSE) {

  if(is.null(marginMatrix)) {

    if(is.null(popTotal)){
      stop("No margin or population total specified for dataMen.")
    }

    writeLines("Calibration only made on population totals for dataMen")
    matrixCal = rep(1,nrow(data))
    formattedMargins = c(popTotal)

  } else {

    # Creation of the elements
    calmarMatrix = marginMatrix[,2:ncol(marginMatrix)]
    # Transform calmarMatrix to numeric matrix to avoid problems in formatMargins
    if(!is.vector(calmarMatrix)) {
      calmarMatrix = matrix(as.numeric(calmarMatrix), nrow=nrow(calmarMatrix), ncol=ncol(calmarMatrix), byrow=F)
    } else {
      calmarMatrix = as.numeric(calmarMatrix)
    }

    popVector <- TRUE
    if(is.null(popTotal)) {
      popVector <- FALSE
    }

    matrixCal = createCalibrationMatrix(marginMatrix,data, popVector)

    formattedMargins = formatMargins(calmarMatrix, matrixCal, popTotal, pct)

  }

  return(list(formattedMargins, matrixCal))

}

## TODO : documentation about integrated calibration
integratedCalibration <- function(dataMen, marginMatrixMen, popMen = NULL,
                                    dataInd, marginMatrixInd, popInd = NULL,
                                    identMen="IDENT_LOG", identInd="IDENT_IND",
                                    colWeights = "POIDS", colCalibratedWeights="POIDS_CALES", method="linear",
                                    maxIter=2500, description=FALSE, bounds, scale=TRUE, check=TRUE) {

  ## TODO : check that idents given are present in tables

  ## Merge right, by identMen
  dataSimultaneous <- merge(dataMen, dataInd, by=identMen, all.y=TRUE)

  #### Count number of ind margin variables by men unit
  # For continuuous variables
  quantiIndMargins <- marginMatrixInd[marginMatrixInd[,2]=="0",]
  if(is.null(ncol(quantiIndMargins))) {
    quantiIndMargins <- quantiIndMargins[1]
  } else {
    quantiIndMargins <- quantiIndMargins[,1]
  }
  print(quantiIndMargins) # debugging
  # For discrete variables
  qualiIndMargins <- marginMatrixInd[marginMatrixInd[,2]!="0",]
  if(is.null(ncol(qualiIndMargins))) {
    qualiIndMargins <- qualiIndMargins[1]
  } else {
    qualiIndMargins <- qualiIndMargins[,1]
  }
  print(qualiIndMargins) # debugging

  ## TODO : quanti variables to dummies
  ## TODO : aggregate individual variables by identMen

  ### Format margins in one margin table
  # Margins for dataMen : TODO -> remove
#   formattedMarginsMen <- createFormattedMargins(dataMen, marginMatrixMen, popMen)[[1]]
#   matrixCalMen <- createFormattedMargins(dataMen, marginMatrixMen, popMen)[[2]]
#   # Margins for dataInd
#   formattedMarginsInd <- createFormattedMargins(dataInd, marginMatrixInd, popInd)[[1]]
#   matrixCalInd <- createFormattedMargins(dataInd, marginMatrixInd, popInd)[[2]]
#
  ## Add individual margins to marginMatrixMen as quantitative variables
  ## TODO


  #### Calibration

  return(dataSimultaneous)
}

##
