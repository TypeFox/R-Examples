## Function to convert a crosstab dataset to long format (database).
##
## Parameters:
## . data      : the crosstab dataset
## . x         : the independent variable to be replicated
## . include   : a character vector containing the names of the columns that should be included
## . exclude   : a character vector containing the names of the columns that should be excluded
## . replicate : a character vector containing the names of columns that have to be replicated for every column name (e.g. variable name)
## . error     : boolean indicating whether the final dataset in long format should contain an extra column for error values (standard deviation, standard error, ...)


cross2long <- function(data, 
                       x,
                       select=NULL,
                       replicate=NULL,
                       error = FALSE,
                       na.rm = FALSE
                      )
{
  # Exception/error handling
  
  varlist <- seq_along(data)
  names(varlist) <- names(data)

  if(!is.null(substitute(x))) {
    indep <- varlist[names(varlist) %in% as.character(substitute(x))]
    indepData <- as.data.frame(data[, indep])
    names(indepData) <- names(varlist)[indep]
  }
  else stop("no independent variable defined")

  reps <- NULL
  if(!is.null(substitute(replicate))) { 
    reps <- varlist[names(varlist) %in% as.character(substitute(replicate))]
    reps <- reps[reps != indep]
    replicatedData <- as.data.frame(data[, reps])
    names(replicatedData) <- names(varlist)[reps]
  }

  if(!is.null(substitute(select))) vars <- varlist[names(varlist) %in% as.character(substitute(select))]
  else vars <- varlist

  vars <- vars[!(vars %in% c(indep,reps))]
  
  selectedData <- as.data.frame(data[, vars])
  names(selectedData) <- names(vars)

  # gathering of data
  nVars <- ncol(selectedData)
  stackedData <- stack(selectedData)
  nObs   <- nrow(stackedData)

  indepData <- as.data.frame(lapply(indepData, FUN=function(x) rep(x,times=nVars)))

  if(!is.null(substitute(replicate))) {
    replication <- as.data.frame(lapply(as.list(replicatedData), FUN=function(x) rep(x,times=nVars)))
    # Construction of appropriate data structure
    if(error)
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indepData[[1]],
                        y    = stackedData$values,
                        err  = rep(1,nObs),
                        as.data.frame(lapply(replication,function(x) factor(x,exclude=NULL)))
                       )
    else
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indepData[[1]],
                        y    = stackedData$values,
                        as.data.frame(lapply(replication,function(x) factor(x,exclude=NULL)))
                       )
  }
  else {
    if(error)
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indepData[[1]],
                        y    = stackedData$values,
                        err  = rep(1,nObs)
                       )
    else
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indepData[[1]],
                        y    = stackedData$values
                       )

  }

  if(na.rm) out <- droplevels(out[complete.cases(out),])

  return(out)
}
