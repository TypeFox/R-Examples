#' Check data for bind function.
#' 
#' check if the data is available for \code{rbind()} or \code{cbind()}
#' 
#' @param data A list containing different sublists ready to be processed by \code{do.call('rbind')} 
#' or \code{do.call('cbind')}
#' @param bind A string showing which bind you are going to use can be 'rbind' or 'cbind'
#' @return data can be processed by bind function; data cannot be processed by bind function
#' @examples
#' data <- list(c(1,1,1),c(2,2,2))
#' bind <- 'rbind'
#' checkBind(data,bind)
#' 
#' data(testdl)
#' \dontrun{
#' checkBind(testdl, 'rbind')
#' }
#' # Since the colnames in testdl are not the same, so it cannot be bound.
#' #
#' @export
checkBind <- function(data, bind){
  # data has to be a list of values, and will be used in do.call('rbind')
  message ('Check if the data list is available for rbind or cbind... \n')
  if (bind == 'rbind') {
    colNum <- sapply(data, function(x) dim(x)[2])
    colLev <- unique(colNum)
    if (length(colLev) != 1) {
      dif <- colLev[2]
      difNum <- which(colNum == dif)
      stop(sprintf('Different Colomn number in %s th of the input list \n', difNum))
      
    }
    
    # For rbind, colnames has  to be checked as well.
    colNameNum <- lapply(data, function(x) colnames(x))
    sameName <- sapply(1:length(colNameNum), function(x) colNameNum[[x]] == colNameNum[[1]])
    if (any(!is.null(unlist(colNameNum))) & (any(sameName == FALSE) | any(length(unlist(sameName)) == 0))) {
      stop('Data in list have Different colnames, which cannot process rbind. ')
    }
    
    
  }else if (bind =='cbind') {
    rowNum <- sapply(data, function(x) dim(x)[1])
    rowLev <- unique(rowNum)
    if (length(rowLev) != 1) {
      dif <- rowLev[2]
      difNum <- which(rowNum == dif)
      stop(sprintf('Different row number in %s th of the input list \n', rowNum))
      
    }
  }
  message('Data list is OK')
}

# Check if a input file is a hyfo grid file.
checkHyfo <- function(...) {
  datalist <- list(...)
  lapply(datalist, FUN = checkHyfo_core)
  invisible()
}

checkHyfo_core <- function(hyfo) {
  #This is to check if the input is a hyfo list.
  checkWord <- c('Data', 'xyCoords', 'Dates')
  if (any(is.na(match(checkWord, attributes(hyfo)$names)))) {
    stop('Input dataset is incorrect, it should contain "Data", "xyCoords", and "Dates",
check help for details or use loadNCDF to read NetCDF file.

If time series input is needed, and your input is a time series, please put "TS = yourinput".')
  }
}

# This check dim is based on the name of the dimension
checkDimLength <- function(..., dim) {
  datalist <- list(...)
  
  for (x in dim) {
    dimLength <- sapply(datalist, function(y) calcuDim(y, x))
    if (any(is.na(dimLength))) stop('No input dimension name, check your dimension name.')
    if (length(unique(dimLength)) != 1) stop('Input data have different dimemsion length.')
  }
  
  invisible()
}




###########################################################################################
##### For biasFactor class

##### Validity functions

checkBiasFactor <- function(object) {
  errors <- character()
  if (length(object@biasFactor) == 0) {
    msg <- 'biasFactors should not be empty.'
    errors <- c(errors, msg)
  }
  
  if (length(object@method) == 0) {
    msg <- 'method should not be empty.'
    errors <- c(errors, msg)
  }
  
  if (length(object@preci) == 0) {
    msg <- 'preci should not be empty.'  
    errors <- c(errors, msg)
  }
  
  prThreshold <- object@prThreshold
  if (length(prThreshold) != 0) {
    if (prThreshold < 0) {
      msg <- 'prThreshold should be greater than 0.'
      errors <- c(errors, msg)
    }
  }
  
  scaleType <- object@scaleType
  if (length(scaleType) != 0) {
    if (scaleType != 'multi' & scaleType != 'add') {
      msg <- paste('scaleType is ', scaleType, '. Should be "multi" or "add".', sep = '')
      errors <- c(errors, msg)
    }
  }
  
  extrapolate <- object@extrapolate
  if (length(extrapolate) != 0) {
    if (extrapolate != 'no' & extrapolate != 'constant') {
      msg <- paste('extrapolate is ', extrapolate, '. Should be "no" or "constant".', sep = '')
      errors <- c(errors, msg)
    }
  }
  
  if (length(errors) == 0) TRUE else errors
}


checkBiasFactor.hyfo <- function(object) {
  errors <- character()
  length_lonLatDim <- length(object@lonLatDim)
  if (length_lonLatDim != 2) {
    msg <- paste('lonLatDim is length ', length_lonLatDim, '. Should be 2', sep = '')
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}


