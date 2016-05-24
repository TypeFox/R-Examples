#' Graphical User Interface for graphical described multiple comparison
#' procedures
#' 
#' Starts a graphical user interface for the creation/modification of directed
#' weighted graphs and applying graphical described multiple comparison
#' procedures.
#' 
#' See the vignette of this package for further details, since describing a GUI
#' interface is better done with a lot of nice pictures.
#' 
#' The GUI can save result files if asked to, can look for a new version on CRAN
#' (if this behaviour has been approved by the user), will change the random seed in 
#' the R session if this is specified by the user in the options (default: no)
#' and could send bug reports if an error occurs and the user approves it.
#' 
#' @param graph Either a variable name for the graph, given as a character
#' string.  (If it is not a syntactically valid name, \code{\link{make.names}}
#' is called to change it to a valid one.)  Or an object of class
#' \code{\link{graphMCP}}.  If the object is modified (even just by updating
#' the class definition or arranging the nodes) it will be saved in the
#' specified environment (default is the global environment).
#' @param pvalues Numeric value that optionally specifies the p-values.
#' @param grid Positive integer that sets the grid size for easier placement of
#' nodes.  (Therefore grid size 1 allows unrestricted placement and disables
#' the grid.)  The default grid=0 uses the last used grid value or if the GUI
#' is started the first time a value of 50.
#' @param debug Logical. If \code{TRUE} debug output is printed to the R
#' console.
#' @param experimentalFeatures Logical. If \code{TRUE} some unfinished /
#' insufficiently tested experimental features are available in the GUI.
#' @param envir Environment where the object \var{graph} is located and/or it
#' should be saved (default is the global environment).
#' @return The function itself returns NULL. But with the GUI a graph can be
#' created or edited that will be available in R under the specified variable
#' name after saving in the specified environment.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @keywords misc graphs
#' @examples
#' 
#' 
#' \dontrun{
#' graphGUI()
#' pvalues <- c(9.7, 1.5, 0.5, 0.6, 0.4, 0.8, 4)/100 
#' graphGUI(HommelEtAl2007(), pvalues=pvalues)
#' 
#' x <- new.env()
#' assign("graph", BonferroniHolm(3), envir=x)
#' graphGUI("graph", envir=x)}
#' 
#' 
#' @export graphGUI
graphGUI <- function(graph="createdGraph", pvalues=numeric(0), grid=0, debug=FALSE, experimentalFeatures=FALSE, envir=globalenv()) {
	if (!is.character(graph)) {
		if ("graphMCP" %in% class(graph)) {
			newGraphName <- "createdGraph"
			i <- 2
			while(exists(newGraphName, envir=envir)) {
				newGraphName <- paste("createdGraph", i, sep="")
				i <- i + 1
			}
			assign(newGraphName, updateGraphToNewClassDefinition(graph), envir=envir)
			graph <- newGraphName
		} else {
			warning("Please specify the variable name for the graph as character.")
			stack <- sys.calls()
			stack.fun <- Filter(function(.) .[[1]] == as.name("graphGUI"), stack)
			graph <- make.names(deparse(stack.fun[[1]][[2]]))
			warning(paste("We guess you wanted to use graphGUI(\"",graph,"\")",sep=""))
		}
	} else {
		if (exists(graph, envir=envir)) {
			if (any(c("graphMCP", "entangledMCP") %in% class(get(graph, envir=envir)))) {
				assign(graph, updateGraphToNewClassDefinition(get(graph, envir=envir)), envir=envir)
				if (is.null(getXCoordinates(get(graph, envir=envir)))||is.null(getYCoordinates(get(graph, envir=envir)))) {
					assign(graph, placeNodes(get(graph, envir=envir)), envir=envir)
				}
			} else {
				stop(paste("The variable",graph,"already exists and is no 'graphMCP' or 'entangledMCP' object."))
			}
		}
	}
  assign("env", envir, envir=gMCPenv)
	invisible(.jnew("org/af/gMCP/gui/CreateGraphGUI", make.names(graph), pvalues, debug, grid, experimentalFeatures))	
}

#' Graphical User Interface for the creation of correlation matrices
#' 
#' Starts a graphical user interface for the correlation matrices.
#' 
#' 
#' @param n Square root of the dimension of the quadratic \eqn{n\times n}{nxn}-Matrix.
#' @param matrix Variable name of matrix of dimension \eqn{n\times n}{nxn} to start with.
#' @param names Row and column names. (Default will be H1,H2,\ldots,Hn.)
#' @param envir Environment where the object \var{matrix} is located and/or it
#' should be saved (default is the global environment).
#' @return The function itself returns NULL.  But with the dialog a symmetric
#' matrix of dimension \eqn{n\times n}{nxn} can be created or edited that will
#' be available in R under the specified variable name after saving.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @keywords misc graphs
#' @examples
#' 
#' \dontrun{
#' corMatWizard(5) # is equivalent to
#' corMatWizard(matrix=diag(5))
#' corMatWizard(names=c("H1", "H2", "H3", E1", "E2"))
#' C <- cor(matrix(rnorm(100),10), matrix(rnorm(100),10))
#' corMatWizard(matrix="C") # or
#' corMatWizard(matrix=C) 
#' }
#' 
#' @export corMatWizard
corMatWizard <- function(n, matrix, names, envir=globalenv()) {  
  if (missing(n) && missing(matrix) && missing(names)) stop("Please specify matrix or dimension.")
  if (!missing(n) && (is.matrix(n) && length(n)>1)) stop("The parameter 'n' should be a single integer number.")
  if (!missing(matrix)) {
    if (!is.character(matrix) || is.matrix(matrix)) {  	
      stack <- sys.calls()
      stack.fun <- Filter(function(.) .[[1]] == as.name("corMatWizard"), stack)
      mname <- deparse(stack.fun[[1]][[2]])    
    }  else {
      mname <- matrix
      matrix <- get(mname, envir=envir)
    }
  } else {
    mname <- "corMat"
    if (!missing(n)) {
      matrix <- diag(n)
    } else {
      matrix <- diag(length(names))
    }
  }
  n <- dim(matrix)[1]  
  if (missing(names)) names <- paste("H",1:n,sep="")
	invisible(.jnew("org/af/gMCP/gui/dialogs/MatrixCreationDialog", as.character(matrix), mname, names))
}
