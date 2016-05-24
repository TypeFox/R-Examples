#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013
# Use granted under BSD license terms
#
# R TFDEA Package
#
# Utilities File
# General internal use only utility functions
#
#******************************************************************************
#
# TODO
# Check for max lp print parameters
# Add debug verbose variables - option
# Log / export info utility - ALL info including orientation, ID string
# * Want to make sure every run is tagged with enough info to reproduce
# Add peers
#
#library(lpSolveAPI)

#
# Package Global Variables
#
tfdea.version   <- "0.9.6"
tfdea.id        <- paste0("tfdea Version ", tfdea.version)
debug           <- 1

#
# Options Strings
#

# DEA
options.orientation.l <- c("input","output")
options.dual.l        <- c(FALSE, TRUE)
options.slack.l       <- c(TRUE, FALSE)
options.second.l      <- c("none", "min", "max")
options.round.l       <- c(FALSE, TRUE)
options.debug.l       <- c(1:4)

# SDEA
options.cook.l        <- c(FALSE, TRUE)

# TFDEA
options.mode.l        <- c("static", "dynamic")


options.rts.l         <- c("vrs","drs", "crs", "irs")
# rts.rhs               <- c(  1,   1,     0,     1)
# rts.typ               <- c( "=", "<=",   0,    ">=")

#
# Set global epsilon for numerical checks in package
# Backup value, dea_internal resets when run from lp control values
epsilon   <- 0.0003162278


# Function: .checkData()
# * Checks for valid types
# * add names if unnamed vars
#
.checkData <- function(x, name){
  # Is value a dataframe, array or vector?
  if ( !is.data.frame(x) && !is.array(x) && (length(x) < 2))
    stop(name," is not a matrix, array (1 or 2 dimensions) or data.frame", call. = FALSE)

  if(length(dim(x)) > 2)
    stop(name," is greater then 2 dimensions", call. = FALSE)

  # If data.frame - convert to matrix - faster to process
  if (is.data.frame(x))
    x <- data.matrix(x)

  # If vector convert to a 1 x N array
  if (is.vector(x))
    x <- array(x, c(length(x), 1))

  # Check that all numeric
  if(! is.numeric(x))
    stop(name," must be numeric", call. = FALSE)

  # Add DMU names if empty
  if (length(dim(x)) >= 1 && is.null(rownames(x)))
    rownames(x) <- paste0("DMU", 1:nrow(x))

  # Add Col names if empty
  if (length(dim(x)) >= 2 && is.null(colnames(x)))
    colnames(x) <- paste0(name, 1:ncol(x))

  return(x)
}

# Function .CheckDataGood
# Check input & output data and warn about problems for DEA
.checkDataGood <- function(x, y, debug=1){

  status <- TRUE

  # Make sure same number of rows
  if (nrow(x) != nrow(y))
    stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)

  # Check for any zero's
  # Check for NA's
  # Check for no positive Values on input & output
  if (any(x == 0, na.rm=TRUE)){
    if (debug >= 1)
      cat("Warning, data has DMU's with inputs that are zero, this may cause numerical problems\n")
    status <- FALSE
  }

  if (any(y == 0, na.rm=TRUE)){
    if (debug >= 1)
      cat("Warning, data has DMU's with outputs that are zero, this may cause numerical problems\n")
    status <- FALSE
  }

  for (k in (1:nrow(x))){
    if ( all(x[k,] <= 0, na.rm=TRUE) & all(y[k,] <= 0, na.rm=TRUE)){
      if (debug >= 1)
        cat("Warning, DMU # ",k, "has no positive inputs or outputs, this may cause numerical problems\n")
    status <- FALSE
    }
  }

  return(status)
}

# Function: .checkVector()
# Check input data and convert to dataframe
# * Checks for valid types
# * add names if unnamed vars
# * make sure returned as array
#
.checkVector <- function(x, name){

  # Is value a dataframe, array or vector?
  if ( !is.data.frame(x) && !is.array(x) && !is.vector(x))
    stop(name," is not a vector, array (1 dimensions) or data.frame", call. = FALSE)

  if(is.data.frame(x))
    x <- data.matrix(x)

  if(!is.null(dim(x))){
    if (length(dim(x)) > 2)
      stop(name," is greater then 2 dimensions", call. = FALSE)
    if (length(dim(x)) == 2 & dim(x)[2] > 1)
      stop(name,"  must have 2nd dimension = 1", call. = FALSE)

    x <- as.vector(x)
  }

  if (!is.numeric(x))
    stop(name, " must be numeric ", call. = FALSE)

  return(x)
}


#
# Check that index is OK, covert to sparse index if logical
.checkIndex <- function(index, x, name){

  if(is.null(index)){                         # NULL special case
    return(c(1:nrow(x)))
  }

  # Is value a dataframe, array or vector?
  if ( !is.data.frame(x) && !is.array(x) && !is.vector(x))
    stop(name," is not a vector, array (1 dimensions) or data.frame", call. = FALSE)

  if(is.data.frame(x))
    x <- data.matrix(x)

# Todo - fix how 1 dim arrays are handled
#   if(!is.null(dim(x))){
#     if (length(dim(x)) > 2)
#       stop(name," is greater then 2 dimensions", call. = FALSE)
#     if (length(dim(x)) == 2 & dim(x)[2] != 1)
#       stop(name," must have 2nd dimension = 1", call. = FALSE)
#
#     x <- as.vector(x)
#   }

  if (is.logical(index)){                     # logical index - convert to sparse index
    if (length(index) != nrow(x))
      stop("Length ",name,"= ",length(index),"  must equal number of DMU's", call. = FALSE)
    return(which(index))
  }

  if (is.numeric(index)){
    if (any(index < 1) | any(index > nrow(x)))
      stop("Value of ",name, "=", paste(index, sep=",", collapse=","),
           " is < 1 or > number of DMU's", call.=FALSE)

    if(anyDuplicated(index)){
      warning("a value of ", name, " is duplicated", call.=TRUE)
    }
  }
  return(index)
}

# Function: .checkOption()
# Check that value in legal list options
# If option in legal list, return lowercase standard form of option
#
.checkOption <- function(value, name, options.l=c(TRUE, FALSE)){

  if (length(value) != 1){
    stop("illegal vector for ", name, " option must be single value", call. = FALSE)
  }

  if (is.character(options.l[1])){                                  # Check character
    if (is.character(value)){
      tmp.value <- tolower(value)
      i <- charmatch(tmp.value, options.l, nomatch=-1)
      if (i > 0)
        return(options.l[i])
    }
  } else if (is.logical(options.l[1])){                             # Check logical
    if (is.logical(value)){
      return(value)
    }
  } else {
    if (is.numeric(value) && is.finite(value) && (value >= 0 )){    # Numeric option
      return(value)
    }
    options.l <- c("numeric >= 0")
  }

  stop("Illegal value=", value, " for ", name, "; legal values are: ",
       paste(options.l, collapse = ", "),
       call. = FALSE)
}


# Check if value is weakly efficient - Depends on orientation
# Needs to be smart FP compare
isEfficient <- function(eff, orientation){

  if (!is.numeric(eff))
    stop("Illegal value for eff; legal values are: numeric", call. = FALSE)

  orientation <- .checkOption(orientation,   "orientation",  options.orientation.l)

  if(orientation == "input"){
    return (is.finite(eff) & (eff + epsilon  >= 1))
  } else {
    return (is.finite(eff) & (eff - epsilon  <= 1))
  }
}

# Internal Function, check if value is efficient Does not depend on orientation
# Needs to be smart FP compare
isStdEfficient <- function(eff){

  # May get called with NaN value - if NaN, return false.

  return( is.finite(eff) & (eff + epsilon) >= 1)
}



#<New Page>
# Function: .lp_debug_file()
#
# Utility function to assist in debugging a linear model.
# Dumps out linear model to a output file.
# Pass a lp_model to print, a log message to use
#
# Use append=FALSE when you want to start a new file, for a new test run for
# example.
#
# NOTE: lp model print function has limit on how big a model it will print.
#       it quietly fails to print the vars if model is too big.
#
# ToDo:
# add a option to control sig digits of values printed
# add a check for if max # vars is being executed and print will fail quietly
# try using another method to dump model out
#
lpDebugFile <- function (lp_model, message="", file_name="lp_model.out", append=TRUE){
  # Check if valid lp object
  if(! (class(lp_model) == "lpExtPtr") )
    stop("not a lp_object", call. = FALSE)

  sink(file_name, type="output", append)  # Send output to file

  cat("\n")
  cat(paste0(message),"\n")
  .print.lpModel(lp_model)
  cat("\n")

  sink(NULL)                              # Resets output to normal location
}

#
# Hacked lpPrint with no limits on size
#
.print.lpModel <- function (x, ...){

  m <- dim(x)[1]
  n <- dim(x)[2]
  control <- lp.control(x)
  if (n < 1) {
    cat(paste("Model name: ", name.lp(x), "\n", sep = ""))
    return(invisible(x))
  }
#   if (n > 8) {
#     cat(paste("Model name: ", name.lp(x), "\n", "  a linear program with ",
#               n, " decision variables and ", m, " constraints\n",
#               sep = ""))
#     return(invisible(x))
#   }
  ans <- matrix(0, m + 1, n)
  for (j in 1:n) {
    col <- get.column(x, j)
    ans[1 + col$nzrow, j] <- col$column
  }
  type <- get.type(x)
  type[type == "integer"] <- "Int"
  type[type == "real"] <- "Real"
  kind <- get.kind(x)
  kind[kind == "standard"] <- "Std"
  kind[kind == "semi-continuous"] <- "S-C"
  bounds <- get.bounds(x)
  upper <- bounds$upper
  lower <- bounds$lower
  ans <- format(rbind(dimnames(x)[[2]], ans, kind, type, upper,
                      lower), justify = "right")
  sense <- ifelse(control$sense == "minimize", "Minimize",
                  "Maximize")
  lhs <- get.constr.value(x, side = "lhs")
  rhs <- get.constr.value(x, side = "rhs")
  rowNames <- format(c("", sense, dimnames(x)[[1]], "Kind",
                       "Type", "Upper", "Lower"))
  constrs <- format(c("", "", get.constr.type(x), "", "", "",
                      ""), justify = "right")
  rhs <- format(c("", "", as.character(rhs), "", "", "", ""),
                justify = "right")
  print.lhs <- any(!is.infinite(lhs[is.element(get.constr.type(x,
                                                               as.char = FALSE), c(1, 2))]))
  lhs <- format(c("", "", as.character(lhs), "", "", "", ""),
                justify = "right")
  if (print.lhs)
    ans <- cbind(rowNames, lhs, constrs, ans, constrs, rhs)
  else ans <- cbind(rowNames, ans, constrs, rhs)
  ans <- apply(ans, 1, paste, collapse = "  ")
  ans <- paste(ans, collapse = "\n")
  model.name <- paste("Model name: ", name.lp(x), "\n", sep = "")
  ans <- paste(model.name, ans, "\n", sep = "")
  cat(ans)
  invisible(x)

}



#
# Utility to print the error message from lp solve for non-zero status
#
lp_solve_error_msg <- function(status){

  error_msg <- c(
    "the model is sub-optimal",
    "the model is infeasible",
    "the model is unbounded",
    "the model is degenerate",
    "numerical failure encountered",
    "process aborted",
    "timeout",
    "the model was solved by presolve",
    "the branch and bound routine failed",
    "the branch and bound was stopped because of a break-at-first or break-at-value",
    "a feasible branch and bound solution was found",
    "no feasible branch and bound solution was found")

  if ( !is.numeric(status) || status < 1 || status > length(error_msg) ){
    warning("Status code:", status, " must be integer from 1 - 12", call. = TRUE)
    message <- paste("ERROR: Unknown status message: ", status)
  } else {
    message <- error_msg[status]
  }
  return(message)
}


.print_status_msg <- function(status, phase, k, debug=1){

  message.error   <- lp_solve_error_msg(status)
  message.std     <- paste("Solver Phase:", phase, "Status:", status, "DMU k=", k)

  if (status != 0){

    if (status == 2 || status == 3){
      if ( debug >= 2) cat(message.std, " is not in technology set: ", message.error, "\n")
    } else {
      if (debug >= 1) cat(message.std, " solver failed: ", message.error, "\n")
    }
  }

}

#
# Error calculations
#

# MAD - Mean Absolute Deviation
# Todo: Add warning messages for NA's
MAD <- function(x, y){
  mad <- mean(abs(x - y), na.rm=TRUE)
  return(mad)
}




