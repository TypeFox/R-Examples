#
# R source code for preprocessing and scoring 
#

No_Causata_Name <- "No Causata Name" # 

#
# constructors for classes
#
CausataVariable <- function(variableName, values, causata.name=variableName) {
  # initialize a Causata Variable object
  # - variableName: a string containing the variable name
  # - values: a vector of values
  
  variableObject <- list()
  class(variableObject) <- "CausataVariable"
  variableObject$name <- variableName   # variable name
  variableObject$causata.name <- causata.name # Name that Causata uses for this variable
  variableObject$class <- class(values) # class of data: numeric, integer, factor
  variableObject$steps <- c() # list of preprocessing steps executed for this variable
  variableObject$outlierUpperLimit <- NULL # upper limit for outlier removal
  variableObject$outlierLowerLimit <- NULL # lower limit for outlier removal
  variableObject$binRight <- NULL # logical, indicating if the intervals should be closed on the right (and open on the left) or vice versa 
  variableObject$binLimits <- NULL # limit values for continous bins
  variableObject$binValues <- NULL # map continuous values to discrete, e.g. weight of evidence binning
  variableObject$missingReplacement <- NULL # missing values will be replaced with this
  variableObject$unmergedLevels <- NULL # levels that remained after MergeLevels
  variableObject$categorizing.attribute <- NULL # Some input variables are categorized.
  variableObject$category.value <- NULL # the category to use for a categorized variable.
  variableObject$time.domain <- NULL # The variable time domain
  
  # transformation is a function that takes a vector (column) from a data.frame and applies
  # all the transformations to it that have been recorded in this CausataVariable, returning the transformed vector
  variableObject$transformation <- identity
  return(variableObject)
}
is.CausataVariable <- function(this) inherits(this, "CausataVariable")


CausataData <- function(dataframe, dependent.variable=NULL, query=NULL) {
  # initialzie a Causata Data object
  # - dataframe: a dataframe containing independent variables
  # - dependent.variable: a vector of dv values, or a name of a variable in the data frame that contains dependent variable
  # - dvName: the name of the dependent variable
  
  dataObject <- list()
  class(dataObject) <- "CausataData"
  
  if (is.null(dependent.variable)){
    # if dependent variable argument is missing then require that "dependent.variable" is defined in data frame
    dataObject$dvName <- "dependent.variable"
    if (!(dataObject$dvName %in% names(dataframe))){
      stop("No argument provided for dependent.variable, so a column named 'dependent.variable' must be provided in the data frame.")
    }
    # copy the dependent variable
    dv <- dataframe$dependent.variable
  } else if (class(dependent.variable) == "character"){
    # character was provided, test if it's a column in the data frame
    if (!(dependent.variable %in% names(dataframe))){
      stop("A variable named ", dependent.variable, " must be present in the data frame.")
    }
    # store the name
    dataObject$dvName <- dependent.variable
    # copy the data to a variable called dv
    dv <- dataframe[[dataObject$dvName]]
  } else {
    # assume that the dependent variable is an array of values with length matching data frame
    # copy the data to a variable called dv
    if (!(nrow(dataframe) == length(dependent.variable))){
      stop("dependent.variable length ", length(dependent.variable), 
           " must match rows in dataframe ", nrow(dataframe))
    }
    if (is.null(names(dependent.variable))){
      # no name provided, so give it a default name
      dataObject$dvName <- "dependent.variable"
    } else {
      # name provided, copy it
      dataObject$dvName <- names(dependent.variable)
    }
    # copy the values to the data frame
    dataframe[[dataObject$dvName]] <- dependent.variable
    dv <- dependent.variable
  }
  
  # dv values can't be missing
  if (any(is.na(dv))){
    stop("Missing values not allowed in dependent variable ", dataObject$dvName)
  }

  # initialize a variable object for every variable in the data frame
  dataObject$variableList <- list()
  dataObject$skippedVariables <- c()
  for (varname in names(dataframe)) {
    # if the dependent variable is in the data frame then skip it, don't create a causata variable for it
    if (varname == dataObject$dvName){
      next # skip to next variable
    }
    # find causata name corresponding to this column of data frame
    #causata.name <- r.to.causata.name.mapping[ match(varname, names(dataframe)) ]
    causata.name <- RToCausataNames(varname)
    if (causata.name == No_Causata_Name){
      # a variable in the data frame does not have a causata name, e.g. profile_id
      # do not create a variable object for this column
      dataObject$skippedVariables <- c(dataObject$skippedVariables, varname)
      next
    } else {
      dataObject$variableList[[varname]] <- CausataVariable(
        varname, dataframe[[varname]], causata.name)
    }
  }
  
  # if variables were skipped then generate a warning
  if (length(dataObject$skippedVariables) > 0){
    warning("CausataData did not create variables for ", length(dataObject$skippedVariables), 
            " that did not conform to variable naming conventions.")
  }
  
  dataObject$df <- dataframe
  dataObject$query <- query
  return(dataObject)
}


# Converts Causata variable names to R variable names.
# For example:
# CausataToRNames(c("total-spend$All Past", "page-view-count$Last 30 Days"))
# returns list("total-spend$All Past"="total.spend__All.Past", "page-view-count$Last 30 Days"="page.view.count__Last.30.Days")
# 
CausataToRNames <- function(name.vector) {
  safe.names <- make.names(str_replace(as.vector(name.vector), "\\$", "__"), unique=TRUE)
  ToRName <- function(index) {
    l <- list()
    l[name.vector[[index]]] <- safe.names[[index]]
    return(l)
  }
  return(unlist(lapply(1:length(name.vector), ToRName)))
}

# Converts R variable column names to Causata variable names.
# For example:
# RToCausataNames(c("total.spend__All.Past", "page.view.count__Last.30.Days"))
# returns list("total.spend_All.Past"="total-spend$All Past", "page.view.count_Last.30.Days"="page-view-count$Last 30 Days")
# 
RToCausataNames <- function(name.vector) {
  DotsToDashes <- function(string) return(str_replace_all(string, "\\.", "-"))
  DotsToSpaces <- function(string) return(str_replace_all(string, "\\.", " "))
  ToCausataName <- function(name) {
    parts <- str_split(name, "__")[[1]] # split at double underscore for data exported from SQL
    parts.csv <- str_split(name, "_")[[1]] # split at single underscore for data exported from csv
    if (length(parts) == 2) {
    #stopifnot(length(parts) == 2)
      return(paste(DotsToDashes(parts[[1]]), DotsToSpaces(parts[[2]]), sep="$", collapse="$"))
    } else if (length(parts.csv)==3 & length(parts.csv[2]>0)) {
      # this variable fits a pattern where there are two underscores with text in between, 
      # assume this is a csv variable, return the name unaltered
      return(name)
    } else {
      # the name doesn't fit pattern name__interval, so set a default non-match value
      warning("The variable named ", name, " does not match Causata naming conventions.")
      return(No_Causata_Name)
    }
  }
  # return a vector, not a list, so use as.vector on results of sapply
  v <- unlist(lapply(name.vector, ToCausataName))
  names(v) <- name.vector
  return(v)
}

#
# define generic functions
#

CleanNaFromContinuous <- function(x, ...) {
  # generic function to replace missing values
  UseMethod("CleanNaFromContinuous")
}

CleanNaFromFactor <- function(x, ...){
  # generic function to replace missing values
  UseMethod("CleanNaFromFactor")
}

Discretize <- function(this, ...){
  # generic function to replace missing values
  UseMethod("Discretize")
}

GetVariable <- function(this, ...) {
  UseMethod("GetVariable")
}

GetTransforms <- function(this, ...) {
  UseMethod("GetTransforms")
}

GetQuery <- function(this, ...) {
  UseMethod("GetQuery")
}


ReplaceOutliers.CausataData <- function(this, variableName, lowerLimit=NULL, upperLimit=NULL, ...){
  # Given a causata data object, replaces outliers with specified values.
  stopifnot(variableName %in% names(this$variableList))
  # extract the causata variable
  causataVariable <- this$variableList[[variableName]]
  # store this step in the list of steps
  causataVariable$steps <- c(causataVariable$steps, "ReplaceOutliers")
  # replace missing values
  this$df[[variableName]] <- ReplaceOutliers( this$df[[variableName]], lowerLimit=lowerLimit, upperLimit=upperLimit )
  # update the variable object
  causataVariable$outlierUpperLimit <- upperLimit
  causataVariable$outlierLowerLimit <- lowerLimit
  this$variableList[[variableName]] <- causataVariable
  # return the data object
  return(this)
}


CleanNaFromFactor.CausataVariable <- function(x, dataObject, replacement, ...) {
  # given a causata variable object, replaces missing values with BLANK
  varname <- x$name
  # store this step in the list of steps
  x$steps <- c(x$steps, "CleanNaFromFactor")
  x$missingReplacement <- replacement
  # replace missing values
  dataObject$df[[varname]] <- CleanNaFromFactor( dataObject$df[[varname]], replacement=replacement )
  # update the variable object
  dataObject$variableList[[varname]] <- x
  # return the data object
  return(dataObject)
}


CleanNaFromFactor.CausataData <- function(x, variableName=NULL, replacement="BLANK", ...) {
  # allow users to pass in a bare, unquoted variable name, this will convert to a string
  if (!is.null(variableName)){
    # a single variable name was provided, confirm that it's in the causataData object
    if (!(variableName %in% names(x$variableList))){
      # the parsed string doesn't match, throw an error
      stop("The variable ", variableName, " was not found in the causataData input.")
    }
    variableList <- variableName
  } else {
    # no variable name was provided, so use all of them
    variableList <- names(x$variableList)
  }

  for (varname in variableList) {
    # test if this is a factor
    if (class(x$df[[varname]]) == "factor") {
      x <- CleanNaFromFactor.CausataVariable( x$variableList[[varname]], x, replacement=replacement )
    }
  }
  # return the data object
  return(x)
}


CleanNaFromContinuous.CausataVariable <- function(x, dataObject, method="median", ...) {
  # given a causata variable object, replaces missing values within
  
  # get the variable name
  varname <- x$name
  # store this step in the list of steps
  x$steps <- c(x$steps, "CleanNaFromContinuous")
  # replace missing values
  outList <- CleanNaFromContinuous( dataObject$df[[varname]], method=method, return.replacement=TRUE )
  # update the data object and store the value used to replace missing values
  dataObject$df[[varname]] <- outList[[1]]
  # check if the missing replacement value was already set, which may have happened in Discretize. If it was set then do not reset here.
  if (is.null(x$missingReplacement)){
    x$missingReplacement <- outList[[2]]
  }
  # update the variable object
  dataObject$variableList[[varname]] <- x
  # return the data object
  return(dataObject)
}


CleanNaFromContinuous.CausataData <- function(x, variableName=NULL, method="median", ...) {
  # allow users to pass in a bare, unquoted variable name, this will convert to a string
  if (!is.null(variableName)){
    # a single variable name was provided, confirm that it's in the causataData object
    if (!(variableName %in% names(x$variableList))){
      # the parsed string doesn't match, throw an error
      stop("The variable ", variableName, " was not found in the causataData input.")
    }
    variableList <- variableName
  } else {
    # no variable name was provided, so use all of them
    variableList <- names(x$variableList)
  }
  
  for (varname in variableList) {
    # test if this is a numeric or integer variable
    varClass <- class(x$df[[varname]])
    if (any(varClass=="numeric") || any(varClass=="POSIXct") || any(varClass=="integer")) {
      x <- CleanNaFromContinuous.CausataVariable(x$variableList[[varname]], x, method=method)
    }
  }
  return(x)
}


Discretize.CausataData <- function(this, variableName, breaks, discrete.values, verbose=FALSE, ...) {
  # Given a causata data object, replaces outliers with specified values.
  stopifnot(variableName %in% names(this$variableList))
  # copy variable from list
  varObject <- this$variableList[[variableName]]
  # require that this variable is numeric
  stopifnot(varObject$class %in% c("numeric","integer"))
  # ensure that missing values were not replaced before this step, missing values should not be replaced before Discretize
  if ("CleanNaFromContinuous" %in% varObject$steps) {
    stop("The CleanNaFromContinuous step must not be executed before the Discretize step for variable ", variableName)
  }
  # ensure that outliers were replaced before this step
  if (!("ReplaceOutliers" %in% varObject$steps)){
    stop("The ReplaceOutliers step must be run before the Discretize step for variable ", variableName)
  }
  # require that the outlier replacement upper limit is <= the highest bin value
  if (is.null(varObject$outlierUpperLimit) || varObject$outlierUpperLimit > breaks[length(breaks)]){
    stop("The outlier limit is not set or greater than the last breaks value in ", variableName)
  }
  if (is.null(varObject$outlierLowerLimit)){
    stop("The outlier lower limit must be set before the Discretize step for variable ",variableName)
  }
  #
  # check the length of inputs and for missing values
  #
  if (any(is.na(this$df[[variableName]]))){
    #
    # missing values are present in numeric variable, must have N+1 breaks and N values
    #
    if (length(breaks) != length(discrete.values)){
      stop("Discretize requires number of breaks to match number of discrete values when missing values are present for variable ", variableName,
           ", found ", length(breaks), " breaks and ", length(discrete.values), " discrete values.")
    }
    # create an artificial bin for missing values, make it above the max value in the data
    newbreak <- 0.5 * (max(breaks) - min(breaks)) + max(breaks)
    breaks <- c(breaks, newbreak)
    # map missing values to the new last break, make it slightly smaller than last artificial bin boundary
    varObject$missingReplacement <- newbreak * 0.99
    this$df[[variableName]] <- CleanNaFromContinuous(this$df[[variableName]], replacement.value = varObject$missingReplacement)
    # update the variable object
    varObject$binLimits <- breaks
    varObject$binValues <- discrete.values
    # set a label to be appended to the last level name for verbose output
    last.level.label.append <- " (Missing)" # append missing label since last level contains missing values
  } else {
    #
    # no missing values present, must have N+1 breaks and N discrete values
    #
    if (length(breaks) != (length(discrete.values)+1)){
      stop("Discretize requires number of breaks N+1 and N discrete values for variable ", variableName, 
           ", found ", length(breaks), " breaks and ", length(discrete.values), " discrete values.")
    }
    # no missing values found, copy over breaks and woe levels as normal
    # update the variable object
    varObject$binLimits <- breaks
    varObject$binValues <- discrete.values
    # set a label to be appended to the last level name for verbose output
    last.level.label.append <- "" # append blank since the last level does not contain missing values
  }
  
  # store this step in the list of steps
  varObject$steps <- c(varObject$steps, "Discretize")
  # set bins to right, closed on right and open on left, default for cut command
  varObject$binRight  <- TRUE
  
  # use cut to discretize the variable
  fiv <- cut(this$df[[variableName]], breaks=breaks, include.lowest=TRUE)
  # replace continuous values with discrete values in the causataData object
  discrete.variable <- rep(0, nrow(this$df))
  level.vec <- levels(fiv)
  if (verbose) {cat("Mapping values to", length(discrete.values), "discrete bins, n is number of values in bins:\n")}
  # loop for each level in the discretized data
  for (i in 1:length(discrete.values)){
    idx <- level.vec[i] == fiv
    if (sum(idx)>0){discrete.variable[idx] <- discrete.values[i]}
    # print message if verbose flag set, add label to level with missing values if present
    if (verbose) {cat("  ", level.vec[i], "->", discrete.values[i], "  n =", sum(idx), 
      if(i==length(discrete.values)) last.level.label.append else "", "\n")}
  }
  this$df[[variableName]] <- discrete.variable
  # update the causata variable and return the causata data object
  this$variableList[[variableName]] <- varObject
  return(this)
}


MergeLevels.CausataVariable <- function(this, causataData, max.levels, ...) {
  # Given a Causata Variable object for a factor, merges the smallest levels 
  # in the factor so that the total number of levels does not exceed maxlevels.
  
  # get the variable name
  varname <- this$name
  # store this step in the list of steps
  this$steps <- c(this$steps, "MergeLevels")
  # Store the levels before merging
  levels1 <- levels(causataData$df[[varname]])
  # merge levels and update the data object
  causataData$df[[varname]] <- MergeLevels(causataData$df[[varname]], max.levels)
  # store a vector of the levels that were not merged into "Other" (all other values become "Other")
  this$unmergedLevels <- levels(causataData$df[[varname]])
  # update the variable object stored in the data object
  causataData$variableList[[varname]] <- this
  # return the data object
  return(causataData)
}


MergeLevels.CausataData <- function(this, variableName=NULL, max.levels, other.name="Other", verbose=FALSE, ...) {
  # given a Causata Data object, merges levels in all of the factors within where
  # the number of levels exceeds maxlevels
  
  if (!is.null(variableName)){
    # a single variable name was provided, confirm that it's in the causataData object
    stopifnot(variableName %in% names(this$variableList))
    variableList <- variableName
  } else {
    # no variable name was provided, so use all of them
    variableList <- names(this$variableList)
  }
  
  if (verbose) {cat("\nMerging levels in factors:\n")}
  for (varname in variableList) {
    # test if this is a factor
    if (class(this$df[[varname]]) == "factor") {
      # it is a factor, check the number of levels
      numLevels <- length(levels(this$df[[varname]]))
      # check if the number of levels exceeds a limit
      if (numLevels > max.levels){
        # this factor exceeds the limit, merge levels
        if (verbose) {cat("   ", varname, "\n")}
      }
      # run mergelevels regardless of whether the threshold number of levels is exceeded.  this way
      # we record which levels were present
      this <- MergeLevels.CausataVariable( this$variableList[[varname]], this, max.levels )
    }
  }
  # return the data object
  return(this)
}


GetVariable.CausataData <- function(this, r.name=NULL, ...) {
  for (variable in this$variableList) {
    if (variable$name == r.name) {
      # exact match found, return variable
      return(variable)
    } else {
      # look for a partial match from a dummy variable produced by model.matrix
      # e.g. var.name__APlevel, where we should match var.name__AP
      idx <- grep( paste("^", variable$name, sep=""), r.name )
      if (length(idx) == 1){
        # match found
        return(variable)
      }
    }
  }
  # if no match was found then throw an error
  stop("No matching variable name found for ", r.name)
}

DiscretizeTransform <- function(this) {
  column.function <- function(column) {
    # if missing values are present then replace them first, before mapping to discrete values
    if (!is.null(this$missingReplacement)){
      # missing replacement value is set, so replace missing values
      column[is.na(column)] <- this$missingReplacement
    }
    this$binValues[ cut(column, breaks=this$binLimits, include.lowest=TRUE, labels=FALSE, right=this$binRight) ]
  }
  ColumnarTransformation(this, column.function)
}

ReplaceOutliersTransform <- function(this) {
  column.function <- function(column) {
    ReplaceOutliers(column, lowerLimit=this$outlierLowerLimit, upperLimit=this$outlierUpperLimit)
  }
  ColumnarTransformation(this, column.function)
}

CleanNaFromContinuousTransform <- function(this) {
  column.function <- function(column) {
    column[is.na(column)] <- this$missingReplacement
    return(column)
  }
  ColumnarTransformation(this, column.function)
}

CleanNaFromFactorTransform <- function(this) {
  ColumnarTransformation(this, function(column) CleanNaFromFactor(column, this$missingReplacement))
}

MergeLevelsTransform <- function(this) {
  stopifnot(is.CausataVariable(this))
  
  unit.of.work <- function(value) {
    if (is.na(value)) {
      NA
    } else if (value %in% this$unmergedLevels) {
      value
    } else {
      "Other"
    }
  }
  
  return(VectorizedFactorTransformation(this, unit.of.work))
}

# Takes a CausataVariable and a function that operates on a single column vector
# returns a function that operates on a data.frame, and applies the given function to the column vector
# that corresponds to the given CausataVariable
#
ColumnarTransformation <- function(this, column.function) {
  function(df) {
    #df[,this$name] <- column.function(df[,this$name])
    # using set from data.table as it is much more efficient with memory
    colnumber <- which(this$name == names(df))
    if (length(colnumber)!=1){
      # require that one column matches
      stop("One column match required but found ", length(colnumber) , " for column ", this$name)
    }
    set(df, j=colnumber, value=column.function(df[, this$name]))
    return(df)
  }
}

# Takes a CausataVariable and a function that operates on a single value in a factor,
# returns a function that operates on a data.frame, and applies the given function to the factor
# that corresponds to the given CausataVariable
#
VectorizedFactorTransformation <- function(this, work.unit, desired.levels=NULL) {
  function(df) {
    col <- df[,this$name]
    stopifnot(is.factor(col))
    new.col <- as.character(col)
    
    for (i in 1:length(new.col)) {
      new.col[i] <- work.unit(new.col[i])
    }
    new.col <- if (length(desired.levels)) {
      factor(new.col, levels=desired.levels)
    } else {
      factor(new.col)
    }
    
    result <- df
    result[this$name] <- new.col
    return(result)
  }
}

# Returns a function that applies the given transforms (functions) in sequence
# Could possibly use Reduce (like a foldleft function)
#
Sequence <- function(functions) {
  function(df) {
    for (f in functions) {
      df <- f(df)
    }
    return(df)
  }
}

GetTransforms.CausataData <- function(this, ...) {
  # iterate over steps, getting a function for each, and the return a function that applies
  # each of those transforms in sequence.
  #
  transforms <- NULL
  for (variable in this$variableList) {
    transforms <- c(transforms, GetTransformFunctionForVariable(variable))
  }
  return(Sequence(transforms))
}

GetQuery.CausataData <- function(this, ...) {
  this$query
}

GetTransformFunctionForVariable <- function(this) {
  stopifnot(is.CausataVariable(this))
  
  transforms <- NULL
  for (step in this$steps) {
    transforms <- c(transforms, switch (step,
          "Discretize"=DiscretizeTransform(this),
          "ReplaceOutliers"=ReplaceOutliersTransform(this),
          "CleanNaFromFactor"=CleanNaFromFactorTransform(this),
          "CleanNaFromContinuous"=CleanNaFromContinuousTransform(this),
          "MergeLevels"=MergeLevelsTransform(this),
          stop("Unknown transformation")
    ))
  }
  return(Sequence(transforms))
}
