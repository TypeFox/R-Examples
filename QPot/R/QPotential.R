#' Computes the quasi-potential for a system of stochastic differential equations using the upwind ordered method. 
#' 
#' 
#' @param x.rhs a string containing the right hand side of the equation for x.
#' @param x.start the starting value of x, usually the x value of the current equilibrium.
#' @param x.bound the x boundaries denoted as c(minimum, maximum).
#' @param x.num.steps the number of steps between the minimum and maximum x value defined in x range.
#' @param y.rhs a string containing the right hand side of the equation for y.
#' @param y.start the starting value of y, usually the y value of the current equilibrium.
#' @param y.bound the y boundaries denoted as c(minimum, maximum).
#' @param y.num.steps the number of steps between the minimum and maximum y value defined in y range.
#' @param save.to.R boolean to output the matrix of results for the upwind-ordered method to the current R session.  The default is to write the matrix to the R session.  save.to.R=FALSE prevents the output from being written to the R session.
#' @param save.to.HD boolean to write the matrix of results for the upwind-ordered method to the hard drive in a file named filename.  Default is FALSE.
#' @param filename string for the name of the file saved to the hard drive.  If save.to.HD=TRUE and filename is left blank, output file saved as defaultname-xX.STARTyY.START.txt, where X.START and Y.START are values in x.start and y.start, respectively.  Matrix stored as a tab-delimited file.
#' @param bounce by default, the upwind-ordered method stops when the boundaries are reached (x.bound and y.bound).  The bounce parameter allows the default action (bounce = 'd'), only positive values to be tested (bounce = 'p'), or reflection near the boundaries (bounce = 'b').
#' @param bounce.edge if bounce = 'b', then to prevent the upwind-ordered method from reaching the boundaries, temporary boundaries are created inside the boundaries defined by x.bound and y.bound.  The boundary edge is bounce.edge of the total range.  Default is 0.01
#' @param k.x integer anisotropic factor for x.  See journal article.  Default is 20.
#' @param k.y integer anisotropic factor for y.  See journal article.  Default is 20.
#' @param verboseC flag (default = TRUE) for printing out useful-for-everyone information from C code implementing the upwind-ordered method (quasipotential.C).
#' @param verboseR NOT IMPLEMENTED: Flag (default = FALSE) for printing out information in QPotential R wrapper.
#' @param debugC NOT IMPLEMENTED: Flag (default = FALSE) for printing out debugging C code 
#' @return filetoHD if save.to.HD enabled, then saves a file in the current directory as either filename or as defaultname-xXSTARTyYSTART.txt
#' @return filetoR if save.to.R enabled, then the function QPotential returns a matrix containing  the upwind-ordered results to be used for plotting.  Requires a variable to catch the returned matrix, i.e. storage <- QPotential(parameters...)
#'
#' @examples
#' # First, System of equations
#' 	equationx <- "1.54*x*(1.0-(x/10.14)) - (y*x*x)/(1.0+x*x)"
#' 	equationy <- "((0.476*x*x*y)/(1+x*x)) - 0.112590*y*y"
#' 
#' # Second, shared parameters for each quasi-potential run
#' 	xbounds <- c(-0.5, 8.0)
#' 	ybounds <- c(-0.5, 8.0)
#' 	xstepnumber <- 200
#' 	ystepnumber <- 200
#' 
#' # Third, a local quasi-potential run
#' 	xinit1 <- 1.40491
#' 	yinit1 <- 2.80808
#' 	storage.eq1 <- QPotential(x.rhs = equationx, x.start = xinit1, 
#'		x.bound = xbounds, x.num.steps = xstepnumber, y.rhs = equationy, 
#'		y.start = yinit1, y.bound = ybounds, y.num.steps = ystepnumber)
#' # Visualize the quasi-potential
#'	QPContour(storage.eq1, dens = c(xstepnumber, ystepnumber), 
#'		x.bound = xbounds, y.bound = ybounds, c.parm = 5) 

QPotential <- function (x.rhs = 'NULL', x.start = 'NULL', x.bound = 'NULL', x.num.steps = 'NULL', y.rhs = 'NULL', y.start = 'NULL', y.bound = 'NULL', y.num.steps = 'NULL', filename = 'NULL', save.to.R = TRUE, save.to.HD = FALSE, bounce = 'd', bounce.edge = 0.01, verboseR = FALSE, verboseC = TRUE, debugC = FALSE, k.x = 20, k.y = 20)
{
# ----------------------------------------------------------------------
# Break apart function parameters into things that C code will use
# ----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Make sure that numeric values are numeric
#-----------------------------------------------------------------------
if ((x.start != 'NULL') && !isTRUE(is.numeric(x.start))) {stop("x.start is not a number")}
if ((x.bound != 'NULL') && !isTRUE(is.numeric(x.bound))) {stop("x.bound is not a number")}
if ((x.num.steps != 'NULL') && !isTRUE(is.numeric(x.num.steps))) {stop("x.num.steps is not a number")}
if ((y.start != 'NULL') && !isTRUE(is.numeric(y.start))) {stop("y.start is not a number")}
if ((y.bound != 'NULL') && !isTRUE(is.numeric(y.bound))) {stop("y.bound is not a number")}
if ((y.num.steps != 'NULL') && !isTRUE(is.numeric(y.num.steps))) {stop("y.num.steps is not a number")}
if ((bounce.edge != 0.01) && !isTRUE(is.numeric(bounce.edge))) {stop("bounce.edge is not a number")}

# ----------------------------------------------------------------------
# check if any component of x is missing
# ----------------------------------------------------------------------
if (verboseR) {print("check if any component of x is missing")}
if (x.rhs == 'NULL') {stop("No equation defined for x. Define x.rhs.")} else {equationx <- x.rhs}
if (grepl(pattern="=", x = x.rhs)) {stop("Equals sign (=) found in x.rhs parameter.  Please give only right hand side of the equation.")}
if (x.start == 'NULL') {stop("No starting value for x. Define x.start.")} else {startxval <- x.start}
if (x.num.steps == 'NULL') {stop("Need the number of steps in x range. Define x.num.steps.")} else {numofstepsx <- x.num.steps}
lengthequationx	<- nchar(x.rhs)
if (length(x.bound) < 2) stop("Not enough values for x range in parameter x.bound.")
if (length(x.bound) > 2) stop("Too many values for x range in parameter x.bound.")
lowerboundsx	<- min(x.bound) #x.bound[1]
upperboundsx	<- max(x.bound) #x.bound[2]

# ----------------------------------------------------------------------
# check if any component of y is missing
# ----------------------------------------------------------------------
if (verboseR) {print("check if any component of y is missing")}
if (y.rhs == 'NULL') {equationy = '0'} else {equationy <- y.rhs}
if (grepl(pattern="=", x = y.rhs)) {stop("Equals sign (=) found in y.rhs parameter.  Please give only right hand side of the equation.")}
lengthequationy <- nchar(y.rhs)
if ( (length(y.bound) == 1) && (y.bound == 'NULL') ) {lowerboundsy <- 0; upperboundsy <- 0}
if (y.bound[1] != 'NULL') {
	if (length(y.bound) < 2) stop("Not enough values for y range in variable y.bound.")
	if (length(y.bound) > 2) stop("Too many values for y range in variable y.bound.")
	lowerboundsy <- min(y.bound) #y.bound[1]; 
	upperboundsy <- max(y.bound) #y.bound[2]
} # end of if y.bound != 'NULL 
if (y.bound[1] == 'NULL') stop('No minimum and maximum y values.  Parameter y.bound not defined')
#The numofstepsy cannot equal 1, because hy=(LY2-LY1)/(NY-1), where NY is numofsteps
if (y.num.steps == 'NULL') {numofstepsy <- 2} else {numofstepsy <- y.num.steps}
if (y.start == 'NULL') {startyval = 0} else {startyval <- y.start}

# ----------------------------------------------------------------------
# Components for C debugging and information to screen
# ---------------------------------------------------------------------- 
if (debugC == FALSE) {is.c.debugging = 0} else {is.c.debugging = 1}
if (verboseC == FALSE) {is.c.loud = 0} else {is.c.loud = 1}

# ----------------------------------------------------------------------
# Components to enable file saving
# ---------------------------------------------------------------------- 
if (verboseR) {print("components to enable file saving")}
#Save in whatever format the user wants
if (filename != 'NULL') {lengthfilename = nchar(filename)}
else 					{lengthfilename = 0}

#default filename has a restriction on the name size
if ( (filename == 'NULL') && (save.to.HD == TRUE) && ( (abs(x.start) > 99999) || (abs(y.start) > 99999) ) ) { stop('Cannot use default filename because program will crash.  Please supply filename')}

if ((save.to.R == TRUE) && (save.to.HD == TRUE)) 		{datasave = 3}
else if ((save.to.R == TRUE) && (save.to.HD == FALSE)) 	{datasave = 2}
else if (isTRUE(save.to.HD))							{datasave = 1}
else													{datasave = 4}
if (save.to.HD == 'testrun')							{datasave = 4}
if (save.to.HD == 'original')							{datasave = 5}
if (verboseR) {print(paste("Variable datasave is: ", datasave, sep = ""))}
#print(paste("datasave = ", datasave,sep=""))

# ----------------------------------------------------------------------
# Determine what C code does at edges of x.bound and y.bound
# bounce (used to make bounce.style), bounce.edge
# ----------------------------------------------------------------------
if (verboseR) {print("Determine what C code does at edges of x.bound and y.bound")}
if (!is.numeric(bounce.edge)) stop('Parameter bounce.edge must be a number')

if (bounce == FALSE) {bounce.style = 'd'}
else if (bounce == 'd') {bounce.style = 'd'}
else if ((bounce == 'p') || (bounce == 'b')) {bounce.style = bounce}
else 	{stop('Parameter bounce must be left blank, (d)efault), (p)ositivevalues, or (b)ounce')}

# ----------------------------------------------------------------------
# Error checking before C code is called
# ----------------------------------------------------------------------

{ 
# check to make sure:
# check that LX1 < xeq < LX2
# check that LY1 < yeq < LY2
if ( (lowerboundsx > startxval) || (upperboundsx < startxval) ) {stop("Starting x value x.start is outside x.bound range")}
if ( (lowerboundsy > startyval) || (upperboundsy < startyval) ) {stop("Starting y value y.start is outside y.bound range")}
if ( lowerboundsx > upperboundsx ) {stop("In x.bound, upper bound is less than lower bound.  x.bound[2] < x.bound[1]")}
if ( lowerboundsy > upperboundsy ) {stop("In y.bound, upper bound is less than lower bounds. y.bound[2] < y.bound[1]")}
 

# ----------------------------------------------------------------------
# warn about memory size as number of steps is very large
# ----------------------------------------------------------------------
if (numofstepsx*numofstepsy > 7000*7000) {
# 8 bytes per numeric, C code and R code each has an array, total number of elements
	tempsizeofarray = 2*8*numofstepsx*numofstepsy
	if (tempsizeofarray > 1e+06) {tempsizeofarraystring <- paste(tempsizeofarray / 1e+06, "MB", sep = " ")}
	if (tempsizeofarray > 1e+09) {tempsizeofarraystring <- paste(tempsizeofarray / 1e+09, "GB", sep = " ")}
	warningstring <- paste("Program expected to use at least", tempsizeofarraystring, "of memory. \n If program crashes, reduce the number of steps and retry.", sep = " ")
	warning(warningstring)
} #end memory size check

} # end of error checking

# ----------------------------------------------------------------------
# Call to upwind ordered method
# ----------------------------------------------------------------------
# Produce multiple versions: 	
#	1) store output only on harddrive 
#	2) store output only in R <DEFAULT>
#	3) store output in R and on harddrive 
#	4) no data saved, testing purposes only
if (datasave == 1) {
	#no R write; HD write
	storage = 0;
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), as.integer(k.x), as.integer(k.y), as.integer(is.c.debugging), as.integer(is.c.loud), PACKAGE="QPot")
#	print(ls())
	return(TRUE)
}
else if (datasave == 2) {
	#DEFAULT
	#R write; no HD write
	storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), as.integer(k.x), as.integer(k.y), as.integer(is.c.debugging), as.integer(is.c.loud), PACKAGE="QPot")
	storage = out2[[1]]
#THIS IS y.num.steps BECAUSE THE C CODE OUTPUTS THE TRANSPOSE
	storage <- matrix(storage, nrow = y.num.steps, byrow = TRUE)
	#1.0e+6 is the INFTY place holder in the C code
	#it means that no QP value was computed
#	tstorage = t(storage)
#	tstorage[tstorage > ((1.0e+6) - 1)] = NA 
#	rm(storage)
#	return(tstorage)
	storage[storage > ((1.0e+6) - 1)] = NA
	return(t(storage))
}
else if (datasave == 3) {
	# R write; HD write
	storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), as.integer(k.x), as.integer(k.y), as.integer(is.c.debugging), as.integer(is.c.loud), PACKAGE="QPot")
	storage = out2[[1]]
	storage <- matrix(storage, nrow = y.num.steps, byrow = TRUE)
	#1.0e+6 is the INFTY place holder in the C code
	#it means that no QP value was computed
#	tstorage = t(storage)
#	tstorage[tstorage > ((1.0e+6) - 1)] = NA 
#	rm(storage)
#	return(tstorage)
	storage[storage > ((1.0e+6) - 1)] = NA
	return(t(storage))
}
else if (datasave == 4) {
	# no R write; no HD write
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), as.integer(k.x), as.integer(k.y), as.integer(is.c.debugging), as.integer(is.c.loud), PACKAGE="QPot")
	return(TRUE)
}
else if (datasave == 5) {
	#uses original c code for file output
	storage = 0;
	print("Note that you have to hand manipulate these data")
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), as.integer(k.x), as.integer(k.y), as.integer(is.c.debugging), as.integer(is.c.loud), PACKAGE="QPot")
}
else {print("datasave is not a possible number.  How did you get here?")}

} # end of QPotential
