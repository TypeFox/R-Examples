#' Inserts parameter values into equations
#'
#' Converts differential equations from string-format (or function-format) with parameters (e.g. "a*x+b) to string-format with parameter values (e.g. 2*x+3).  Specifically, Model2String reads in the equations, searches for the differential equations within the function (if required), and replaces the parameters with numerical values given by the user.  Returns an array of strings containing the differential equations.
#' 
#' @param model contains the differential equations as given to \code{\link{TSTraj}}.  Can either be a string or a function used by the package \code{deSolve} (see third example).
#' @param parms a named vector of paramters and their respective values for the deterministic equations.  If inputing a function and parms is empty, Model2String will return the equation as a string.
#' @param deSolve.form if FALSE (default) interprets model as a string containing the rhs of the equation.  If TRUE, treats model as the function input into the package \code{deSolve}.
#' @param x.lhs.term string containing the left hand side of the first equation to search for if deSolve.form is TRUE. Default is 'dx'.
#' @param y.lhs.term string containing the left hand side of the second equation to search for if deSolve.form is TRUE.  Default is 'dy'.
#' @param supress.print if FALSE (default), suppress output.  TRUE prints out equations converted to strings with parameters replaced with values.
#' @param width.cutoff parameter width.cutoff from \code{deparse} in package:base.  Determines the cutoff in bytes at which line breaking is tried.  Default is 500 with possible range of 2 to 500.
#' @return equations: an array with strings of the equations with substituted parameter values.  If deSolve.form is TRUE, the first is the x equation, the second is the y equation.
#'
#' @examples

#' # First example with the right hand side of an equation
#' test.eqn.x = "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))"
#' model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1)
#' equations.as.strings.x <- Model2String(test.eqn.x, parms = model.parms)
#'
#' # Second example with individual strings with left and right hand sides
#' # Note the use deSolve.form = TRUE to remove the dx and dy
#' test.eqn.x = "dx = (alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))"
#' test.eqn.y = "dy = ((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)"
#' model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1, gamma=0.476, mu=0.112509)
#' equations.as.strings.x <- Model2String(test.eqn.x, parms = model.parms, 
#'  deSolve.form = TRUE, x.lhs.term = 'dx', y.lhs.term = 'dy')
#' equations.as.strings.y <- Model2String(test.eqn.y, parms = model.parms, 
#'  deSolve.form = TRUE, x.lhs.term = 'dx', y.lhs.term = 'dy')
#'
#' # Third example with deSolve-style function call:
#' model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1, gamma=0.476, mu=0.112509)
#' ModelEquations <- function(t, state, parms) {
#' 	with(as.list(c(state, parms)), {
#'	dx <- (alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))
#'	dy <- ((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)
#'	list(c(dx,dy))
#'	})
#' }
#'
#' Model2String(ModelEquations, parms = model.parms, deSolve.form = TRUE, 
#'  x.lhs.term = 'dx', y.lhs.term = 'dy') 


Model2String <- function(model, parms = 'NULL', deSolve.form = FALSE, x.lhs.term = 'dx', y.lhs.term = 'dy', supress.print = FALSE, width.cutoff = 500) {
	if (!supress.print) {
		print("Note: This function is supplied as duct tape.  Long equations, equations spanning multiple lines, equations with strange notation, etc, may not work.  Always check the output.")
	}
	
	if (deSolve.form == TRUE) {

		temp <- deparse(model, width.cutoff = width.cutoff)
	# for some reason, inputing a string causes a " to be added
		temp <- gsub(pattern = '\"', replacement = '', x = temp)
	
	#Dump function into a list of character strings
	#Go through each string and determine if it contains an equation
		foundx = 0	#flag for making sure dx is only found once
		foundy = 0	#flag for making sure dy is only found once
	
	#remove the lhs and return the rhs
		for (i in 1:length(temp)) {
		#when searching, first look for the lhs defining whether the derivative is for x or y
		#once found, look inside the string and use either '<-' or '=' to separate
		# the lhs from the rhs
			if ((foundx == 0) && isTRUE(grepl(pattern = x.lhs.term, x = temp[i]))) {
				foundx = 1
				if (grepl(pattern = '<-', x = temp[i])) {
					location <- regexpr(pattern = '<-', text = temp[i])
					x.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
				}
				else if (grepl(pattern = '=', x = temp[i])) {
					location <- regexpr(pattern = '=', text = temp[i])
					x.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
				} else {stop("Equation does not contain = or <-")}
			}
	
			if ((foundy == 0) && isTRUE(grepl(pattern = y.lhs.term, x = temp[i]))) {
				foundy = 1
				if (grepl(pattern = '<-', x = temp[i])) {
					location <- regexpr(pattern = '<-', text = temp[i])
					y.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
				}
				else if (grepl(pattern = '=', x = temp[i])) {
					location <- regexpr(pattern = '=', text = temp[i])
					y.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
				} else {stop("Equation does not contain = or <-")}
			}
		}
		
		equations = c()
		if (foundx == 1) {equations = x.equation}
		if (foundy == 1) {equations = c(equations, y.equation)}
	} #if (deSolve.form == TRUE)
	else { #simple find and replace
		equations <- model
	}
	
#if parameters are not declared, then we do not have to replace anything
	if (!(parms[1] == 'NULL')){
#replace the parameter names in the equations with their values
		allnames <-names(parms)
		for (i in 1:length(parms)) {
			currname <- allnames[i]
			value = toString(parms[[i]])
			equations <- gsub(pattern = currname, replacement = value, x = equations)
		}
	} #if (!(parms[1] == 'NULL')){
	
	if (!supress.print) {print(paste(equations))}
	return(equations)

}
