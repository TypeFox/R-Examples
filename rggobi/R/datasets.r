# Basic data set manipulation
#
# Subsetting, dimensions, names etc.
# Attempts to mimick R data frames as closely as possible
# ============================================================================

# Print GGobiData
# Print GGobiData
#
# By default printing a GGobiData acts like
# printing an R data.frame - ie. show all the data
#
# @arguments GGobi dataset to retrieve
# @keyword attribute
# @keyword internal
print.GGobiData <- function(x, ...) {
  print(as.data.frame(x), ...)
}

# GGobiData dimensions
# Retrieve the dimension of a GGobiData
#
# @arguments dataset
# @keyword attribute
# @keyword internal
#X g <- ggobi(mtcars)
#X dim(g[1])
dim.GGobiData <- function(x) {
	.GGobiCall("datasetDim", x)
}

# GGobiData rows
# Retrieve the number of row in a GGobiData
#
# @arguments dataset
# @keyword attribute
# @keyword internal
nrow.GGobiData <- function(d) dim(d)[2]

# GGobiData columns
# Retrieve the number of columns in a GGobiData
#
# @arguments dataset
# @keyword attribute
# @keyword internal
ncol.GGobiData <- function(d) dim(d)[1]

# GGobiData column names
# Get column names for a GGobiData
#
# @arguments dataset
# @keyword attribute
# @keyword internal
names.GGobiData <- function(x, ...) {
  .GGobiCall("getVariableNames", FALSE, x)
}

# Set column names
# Set column names for a GGobiData
#
# @arguments GGobiData
# @arguments new names
# @keyword attribute
# @keyword internal
"names<-.GGobiData" <- function(x, value) {
	colnames(x) <- value
	x
}


# Variable index
# Return indices corresponding to variable names
#
# @arguments GGobiData
# @arguments variable names
# @returns numeric vector of (0-based) positions
# @keyword attribute
# @keyword internal
variable_index <- function(x, names) {
	if(length(names) == 0) return(integer(0))
	#if(is.integer(names)) return(as.integer(names - 1))
	if(is.numeric(names)) return(as.integer(names - 1))
	if(is.character(names)) return(as.integer(match(names, names(x)) - 1))

	def <- list(X = numeric(0), Y=numeric(0), Z=numeric(0))
	if (!is.null(names$X)) def$X <- names$X
	if (!is.null(names$Y)) def$Y <- names$Y
	if (!is.null(names$Z)) def$Z <- names$Z

	lapply(def, function(name) variable_index(x, name))
}

# Get row names
# Get row names for a GGobiData
#
# @arguments ggobiDataget
# @arguments new names
# @keyword attribute
# @keyword internal
rownames.GGobiData <- function(x) {
	.GGobiCall("getRowNames", x)
}

# Set row names
# Set row names for a GGobiData
#
# @arguments GGobiData
# @arguments new names
# @keyword attribute
# @keyword internal
#X g <- ggobi(mtcars)
#X df <- g[1]
#X rownames(df)
#X rownames(df) <- tolower(rownames(df))
#X rownames(df)
"rownames<-.GGobiData" <- function(x, value) {
	dims <- dimnames(x)
	stopifnot(length(x) == length(dims[1]))
	dims[1] <- value
	dimnames(x) <- dims
	x
}

# Get dimension names
# Get row and column names for a GGobiData
#
# @arguments ggobiDataget
# @keyword attribute
# @keyword internal
dimnames.GGobiData <- function(x) {
  list(rownames.GGobiData(x), names(x))
}

# Set dim names
# Set dim names for a GGobiData
#
# @arguments GGobiData
# @arguments new names
# @keyword attribute
# @keyword internal
"dimnames<-.GGobiData" <- function(x, value) {
  .GGobiCall("setRowNames", as.character(value[[1]]), as.integer(1:length(value[[1]]) - 1), x)
  .GGobiCall("setVariableNames", as.integer(1:ncol(x) -1 ), as.character(value[[2]]), x)
  x
}

# Summarise GGobiData.
# Summarise a GGobiData with dimensions, mode and variable names.
#
# @arguments GGobiData
# @arguments ignored
# @keyword attribute
summary.GGobiData <- function(object, ...) {
	list(dim = dim(object), mode = mode(object), variables = names(object))
}

# Subsettting
# Subsetting for ggobi datasets
#
# This functions allow one to treat a ggobi dataset as if it were a local
# data.frame.  One can extract and assign elements within the dataset.
#
# This method works by retrieving the entire dataset into
# R, and then subsetting with R.
#
# @arguments ggobi dataset
# @arguments rows
# @arguments cols
# @arguments drop dimensions?
# @value desired subset from data.frame
# @alias [[.GGobiData
# @alias $.GGobiData
# @keyword manip
#X g <- ggobi(mtcars)
#X x <- g$mtcars
#X x[1:5, 1:5]
#X x[[1]]
#X x$cyl
"[.GGobiData" <- function(x, i, j, drop=FALSE) {
	as.data.frame(x)[i, j, drop=drop]
}

"$.GGobiData" <- "[[.GGobiData" <- function(x, i) {
	as.data.frame(x)[[i]]
}


# Conversion methods
# Convert a GGobiData to a regular R data.frame or matrix
#
# @arguments GGobiData
# @keyword manip
# @keyword internal
# @alias as.matrix.GGobiData
"as.data.frame.GGobiData" <- function(x, ...) {
	gd <- .GGobiCall("getData", x)
	df <- as.data.frame(gd)
	rownames(df) <- rownames(x)
	df
}

"as.matrix.GGobiData" <- function(x, ...) {
  as.matrix(as.data.frame(x))
}

# Assignments for ggobi datasets
# This functions allow one to treat a ggobi dataset as if it were a local data.frame.  One can extract and assign elements within the dataset.
#
# This method works by retrieving the entire dataset into
# R, subsetting that copy, and then returning any changes.
#
# @arguments ggobi dataset
# @arguments row indices
# @arguments column indices
# @arguments new values
# @keyword manip
# @keyword internal
# @alias $<-.GGobiData
# @alias [[<-.GGobiData
#X g <- ggobi(mtcars)
#X x <- g["mtcars"]
#X x[1:5, 1:5]
#X x[1:5, 1] <- 1:5
#X x[1:5, 1:5]
"[<-.GGobiData" <- function(x, i, j, value) {
	data <- as.data.frame(x)

	# figure out if any new columns have been added and add them.

  if (missing(i))
    data[,j] <- value
  else if (missing(j))
    data[i,] <- value
  else data[i, j] <- value
	for(var in unique(j)) {
		x[[var]] <- data[[var]]
	}
	x
}

"$<-.GGobiData" <- "[[<-.GGobiData" <- function(x, i, value) {
	df <- as.data.frame(x)
	if (is.null(value)) {
		ggobi_data_remove_variable(x, i)
	} else if (is.null(df[[i]])) {
		ggobi_data_add_variable(x, value, i)
	} else {
		ggobi_data_set_variable(x, value, i)
	}
	x
}


# Set variable values
# Set the variable values for a column in a GGobiData
#
# @arguments GGobiData
# @arguments values of new variable
# @arguments variable name
# @arguments update?
# @keyword internal
ggobi_data_set_variable <- function(x, vals, var, update = TRUE) {
	varId <- variable_index(x, var)
	if(any(is.na(varId))) stop("Invalid variable")

	.GGobiCall("setVariableValues", as.numeric(vals), as.integer(1:length(vals) - 1), varId, as.logical(update), x)
	x
}

# Add variable
# Add variable to a GGobiData
#
# @alias addVariable
# @arguments GGobiData
# @arguments values to add
# @arguments name of column to add
# @keyword internal
ggobi_data_add_variable <- function(x, vals, name, ...) {
	if (!(is.factor(vals) || is.numeric(vals))) stop("Variable must be a factor or numeric")
	if (length(vals) != nrow(x)) stop("Variable must be same length as existing data set")

	levels <- NULL
	values <- NULL

	if(is.factor(vals)) {
		levels <- table(vals)
		values <- sort(unique(as.integer(vals)))
	}

	.GGobiCall("addVariable", as.numeric(vals), as.character(name), levels, values, x)
}

# Remove variable
# Remove variables from a GGobiData object
#
# @keyword internal
ggobi_data_remove_variable <- function(x, var) {
	varId <- variable_index(x, var)
	if(any(is.na(varId))) stop("Invalid variable")

 .GGobiC("removeVariables", varId, x)[[1]]
}


ids <- function(x) UseMethod("ids", x)
"ids<-" <- function(x, value) UseMethod("ids<-", x)

# Row ids
# Retrive row ids from a GGobiData
#
# @alias ids
# @arguments GGobiData
# @keyword manip
# @seealso \code{\link{ids<-}}
ids.GGobiData <- function(x) {
  .GGobiCall("getCaseIds", x)
}

# Set row ids
# Set row ids from a GGobiData
#
# @alias ids<-
# @arguments GGobiData
# @arguments new values
# @keyword manip
# @seealso \code{\link{ids}}
"ids<-.GGobiData" <- function(x, value) {
	.GGobiCall("setIDs", as.character(value), .data=x)
	x
}

## for S4 integration
# setOldClass("GGobiData")
