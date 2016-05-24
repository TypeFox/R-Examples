# Display construction
# Create a new display
# 
# @keyword internal 
display <- function(x, ...) UseMethod("display", x)

# Create a new display
# Create a new display for the GGobiData object.  
#
# This function allows you to create a new display from
# a GGobiData object.  You will need to specify the type of display you 
# want ("Scatterplot Display", "Scatterplot Matrix" and  "Parallel Coordinates
# Display" are the most common), and which variables the plot
# should be intialised with.  Specifying more than two variables only makes
# sense for scatterplot matrices and pcps.
# 
# Many of the plots used in GGobi (eg. the tours and densities plots) are
# special modes of the scatterplot display.  You will need to create a 
# new scatterplot display, change the projection mode to what you want, 
# and then set the variables.  Hopefully this will be improved in a future
# version of rggobi.
#  
# @arguments GGobiData object
# @arguments projection mode to use
# @arguments variables to display, see \code{\link{variables.GGobiDisplay}} for details
# @arguments If TRUE, returns widget for use with RGtk2
# @arguments ignored
# @seealso \code{\link{ggobi_display_types}} for a list of display types
# @keyword dynamic
#X g <- ggobi(mtcars)
#X display(g[1])
#X display(g[1], vars=list(X=4, Y=5))
#X display(g[1], vars=list(X="drat", Y="hp"))
#X display(g[1], "Parallel Coordinates Display")
#X \dontrun{
#X display(g[1], "2D Tour")
#X display(g[1], "2x1D Tour", list(X=c(1,2,3), Y=c(4,5,6)))
#X } 
#X display(g[1], "Scatterplot Matrix")
display.GGobiData <- function(x, pmode="Scatterplot Display", vars=list(X=names(x)), embed=FALSE, ...) 
{
	type <- pmodes()[pmode]
	ivars <- variable_index(x, vars)
  embed <- as.logical(embed)

	d <- .GGobiCall("createDisplay", ggobi_display_make_type(type), ivars$X, x, !embed)
	if (type != pmode) {
		pmode(d) <- pmode
	}
  variables(d) <- vars
	d
}


# Get variables
# Get variable names from a container
#
# @keyword internal 
variables <- function(x) UseMethod("variables", x)

# Set variables
# Set variable names in a container
#
# @keyword internal 
"variables<-" <- function(x, value) UseMethod("variables<-", x)

# Get display variables
# List the variables used in a given display
# 
# There are three types of variables in GGobi displays: 
# X, Y, Z, which correspond to the labels on the toggle buttons
# in GGobi.  Most plots have a constrained set of possible options.
# For example, in tours you can only set X variables, and you must 
# have at least three.  Or in the rotation plot, you need exactly
# one X, Y, and Z variable.   You can figure out what these 
# conditions are by using the toggle buttons in GGobi.
# 
# @arguments GGobiDisplay object
# @keyword dynamic 
# @seealso \code{\link{variables<-.GGobiDisplay}} for examples
variables.GGobiDisplay <- function(x) {
  vars <- .GGobiCall("getDisplayVariables", x, .gobi=ggobi(x))

  names(vars[[1]]) <- colnames(dataset(x))[vars[[1]] + 1]
  split(vars[[1]] + 1, factor(vars[[2]], levels = c("X", "Y", "Z")))
}

# Set display variables 
# Set display variables with a list of x, y, and z component variable indices.
# 
# There are three types of variables in GGobi displays: 
# x, y, z, which correspond to the labels on the toggle buttons
# in GGobi.  Most plots have a constrained set of possible options.
# For example, in tours you can only set x variables, and you must 
# have at least three.  Or in the rotation plot, you need exactly
# one x, y, and z variable. 
#
# Currently, there is no checking done to ensure that you are
# sending a sensible set of variables for the given display type.
# Generally, any invalid choices will be silently ignored.
# 
# @arguments GGobiDisplay object
# @arguments list with X, Y and Z components listing the variable indices to display, either as numeric position or character variable name
# @keyword dynamic 
#X g <- ggobi(mtcars)
#X d <- display(g[1], "Parallel Coordinates Display")
#X variables(d)
#X variables(d) <- list(X=1:8)
#X variables(d) <- list(X=c("mpg", "cyl"))
#X variables(d)
"variables<-.GGobiDisplay" <- function(x, value) {
	stopifnot(is.list(value))
	names(value) <- toupper(names(value))

	d <- dataset(x)
	prev_vars <- variable_index(d, variables(x))
	new_vars <- variable_index(d, value)

  if (any(is.na(new_vars)))
    stop("Variable(s) ", paste(value[is.na(new_vars)], collapse=", "), 
      " not in dataset")
  
  prev_vars <- mapply(setdiff, prev_vars, new_vars, SIMPLIFY=FALSE)

	.GGobiCall("setDisplayVariables", new_vars, prev_vars, x, .gobi=ggobi(x))
	x
}
