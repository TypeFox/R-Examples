# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## utilities for class "BasicVector"

# get character eqivalent
getCharacter <- function(x, names) {
    if(is(x, "character")) x
    else if((is(x, "logical") || is(x, "numeric")) && length(x)) names[x]
    else character()  # other classes
}

# get length of specified selection
getSelectionLength <- function(x) {
    if(is(x, "character")) length(x)
    else if(is(x, "numeric") && all(x >= 0)) length(x[x > 0])
    else NA  # other classes
}

# ---------------------------------------

## utilities for class "NumericMatrix"

# check if data is numeric
checkNumericMatrix <- function(x) {
    if(is(x, "numeric")) TRUE 
    else if(is(x, "matrix")) mode(x) == "numeric"
    else FALSE  # other classes
}

# ---------------------------------------

## utilities for simulation results

# get NA rate for simulation results
convertNArate <- function(x) {
    if(is(x, "numeric")) x
    else if(is(x, "matrix")) 1:nrow(x)
    else numeric()  # other classes
}

# get number of repetitions for additional information
getRepetitions <- function(x) {
    if(is(x, "numeric")) if(length(x)) 1 else 0 
    else if(is(x, "data.frame") || is(x, "matrix")) nrow(x)
    else numeric()  # other classes
}

# ---------------------------------------

## utilities for plot functions

# get formula for plot functions
getFormula <- function(left, right, cond) {
    if(length(cond)) {
        cond <- paste(cond, collapse=" + ")
        as.formula(paste(left, "~", right, "|", cond))
    } else as.formula(paste(left, "~", right))
}

# get data in the correct format for lattice graphics
getLatticeData <- function(x, cond, names) {
    n <- nrow(x)
    tmp <- lapply(names, 
        function(nam) {
            cbind(x[, cond, drop=FALSE], .Name=rep.int(nam, n), .Value=x[, nam])
        })
    do.call(rbind, tmp)
}

# ---------------------------------------


## other utilities

# get argument names of a function
argNames <- function(fun, removeDots = TRUE) {
    nam <- names(formals(fun))
    if(removeDots) nam <- setdiff(nam, "...")
    nam
}

# drop unused levels in case of factor, convert otherwise
getFactor <- function(x) {
    if(is.factor(x)) x[, drop=TRUE]
    else as.factor(x)
}

# get names of real columns of a data.frame 
# (i.e., remove those used internally by simFrame)
getNames <- function(x) {
#    nam <- names(x)
#    nam[substring(nam, 1, 1) != "."]
    setdiff(names(x), c(".weight",".contaminated"))
}

# remove internal variables from the data.frame
removeInternal <- function(x, drop = FALSE) x[, getNames(x), drop=drop]
