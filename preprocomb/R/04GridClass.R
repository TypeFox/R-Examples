#' @include 03PreprocessorClass.R
NULL

reportexitstatus <- function(preprocesseddatasets){

  if(class(preprocesseddatasets)!="list"){stop("Argument 'preprocesseddatasets' to function 'reportexitstatus' must of a list.")}

  variance <- unlist(lapply(preprocesseddatasets, function(x) slot(x, "variance")))
  finite <- unlist(lapply(preprocesseddatasets, function(x) slot(x, "finite")))
  completeobs <- unlist(lapply(preprocesseddatasets, function(x) slot(x, "completeobs")))
  classbalance <- unlist(lapply(preprocesseddatasets, function(x) slot(x, "classbalance")))
  ntopratiotwoplus <- unlist(lapply(preprocesseddatasets, function(x) slot(x, "ntopratiotwoplus")))
  mindimensions <- unlist(lapply(preprocesseddatasets, function(x) slot(x, "mindimensions")))

  check <- all(c(variance, finite, completeobs, classbalance, ntopratiotwoplus, mindimensions))==TRUE

  if (check==TRUE) {exitstatus <- c("Exit status: OK: Stable computation of misclassification errors expected.")}
  if (check==FALSE) {exitstatus <- c("Exit status: Warning: Unstable computation of misclassification errors expected. See: yourgridclassobject@data")}

  return(exitstatus)

}

creategrid <- function(phases){
  if(class(phases)!="list"){stop("Argument 'phases' to function 'creategrid' must of a list.")}
  grid <- expand.grid(lapply(phases, function(x) eval(as.name(x))@preprotransformations))
  colnames(grid) <- unlist(phases)
  return(grid)
}

## DATA FORMATION

initializedataslot <- function(classname, dataobject){

  tryCatch({

  subclassobject <- new(classname)

  if (class(dataobject)=="DataClass") {transformeddata <- transformdata(subclassobject, dataobject)}

  if (class(dataobject)=="data.frame") {transformeddata <- transformdata(subclassobject, initializedataclassobject(dataobject))}

  if (is(dataobject, "PreprocessorClass")==TRUE) {transformeddata <- transformdata(subclassobject, dataobject@data)}

  subclassobject@data <- transformeddata

  return(subclassobject)

  }, error= function(e) return({subclassobject <- new(classname)})
  )

}

executepreprocessing <- function(grid, dataobject){

  if(class(grid)!="data.frame"){stop("Argument 'grid' must of a data frame of a GridClass object.")}
  if(class(dataobject)!="DataClass"){stop("Argument 'dataobject' must of a DataClass object.")}

  out_preprocesseddatasets <- vector(mode="list", nrow(grid))
  firstcolumningrid <- 1

  for (rowingrid in 1:nrow(grid))
  {
    out_preprocesseddatasets[rowingrid] <- initializedataslot(as.character(grid[rowingrid, firstcolumningrid]), dataobject) # first column of grid

    if (ncol(grid) > 1){

      for (columningrid in 2:ncol(grid))
      {
        out_preprocesseddatasets[rowingrid] <- initializedataslot(as.character(grid[rowingrid,columningrid]), out_preprocesseddatasets[[rowingrid]]@data)
      }

    }

  }

  out_preprocesseddatasets <- lapply(out_preprocesseddatasets, function(x) slot(x, "data"))
  out_preprocesseddatasets <- lapply(out_preprocesseddatasets, validatedata)
  print(reportexitstatus(out_preprocesseddatasets))
  return(out_preprocesseddatasets)
}


## GRID

#' GridClass
#'
#' GridClass is a container for preprocessor combinations and the corresponding preprocessed data sets.
#' GridClass is an interface for extending the system.
#' @slot grid (data frame) preprocessor combinations
#' @slot data (list) DataClass objects
#' @slot validation (data frame) validation results
#' @details Extensions can include approximate combinatorial optimization for finding near-best
#' combinations faster.
#' @export

setClass("GridClass", representation(grid="data.frame", data="list", validation="data.frame"))

#' setgrid
#'
#' setgrid takes the preprocessing phases, which contain preprocessors and creates
#' the combinations of them as a grid. It then computes and stores the transformed
#' data sets for each combination. setgrid initializes a GridClass object.
#
#' @param phases (character) vector of phases
#' @param data (data frame)
#' @param diagnostics (logical) run testpreprocessor(), defaults to TRUE
#' @return a GridClass object
#' @examples
#' ## grid <- setgrid(phases=c("outliers", "selection"), data=iris)
#' @details If there are missing value, imputation phase must be set as first phase.
#' @export

setgrid <- function(phases, data, diagnostics=TRUE){

# Validate arguments
if(class(phases)!="character"){stop("Argument 'phases' must be a character vector.")}
if(class(data)!="data.frame"){stop("Argument 'data' must of a data frame.")}

phases <- as.list(phases)
if(!all(lapply(phases, function(x) class(eval(as.name(x))))=="PhaseClass")){
stop("All elements in argument 'phases' must point to PhaseClass objects.")}

# Initialize objects

dataclassobject <- initializedataclassobject(data)

gridclassobject <- new("GridClass")

gridclassobject@grid <- creategrid(phases)

# Validate preprocessors

if (diagnostics==TRUE){
print("Running diagnostics on single preprocessors:")
validation <- testpreprocessors(unique(unlist(gridclassobject@grid)))
}
print("Preprocessing data set by combinations:")
gridclassobject@data <- executepreprocessing(gridclassobject@grid, dataclassobject)

# Validation results on combinations
validationresults <- data.frame(gridclassobject@grid, data.frame(t(data.frame(lapply(gridclassobject@data, extract)))))
row.names(validationresults) <- NULL
gridclassobject@validation <- validationresults

return(gridclassobject)
}




