#' Turn a Character String of the Right-Hand Side of an Equation into an R Function
#'
#' \code{equation2function} converts the right-hand side of an equation that can be written
#' as a single line (like the right-hand side of an object of class \code{formula}) and
#' creates an R function which calls the variables within as arguments and returns values
#' consistent with the parameters of the input equation as written.

#' @details
#' This simple little function is rather hacky but seems to get the job done, for a
#' functionality that doesn't seem to be present elsewhere in R.

#' @param equation The right-hand-side (RHS) of an equation, given as a character string.
#' If not of type character, \code{equation2function} attempts to coerce \code{equation} to
#' type character and fails if it cannot.

#' @param envir The environment the resulting function will be evaluated in.
#' See \code{\link{as.function}}.

#' @param notName A useless string used simply as a placeholder in turning \code{equation}
#' into a function, which should not match any actual variable in \code{equation}. Only
#' supplied as an argument in case any 

#' @return
#' A function, with named blank (i.e. no default value) arguments.

#' @author David W. Bapst

#' @examples
#' # some simple examples
#' foo<-equation2function("x+y")
#' foo
#' foo(x=4,y=0.1)
#'
#' foo<-equation2function("x+2*sqrt(2*y+3)^2")
#' foo
#' foo(x=4,y=0.1)
#'
#' # what about weird long argument names and spaces
#' foo<-equation2function("stegosaur + 0.4 * P")
#' foo
#' foo(stegosaur=5,P=0.3)

#' @name equation2function
#' @rdname equation2function
#' @export
equation2function<-function(equation,envir=parent.frame(),
	notName="XXXXXXXXXXX"){
	#checks
	if(!is.character(equation)){
		equation<-as.character(equation)
		if(!is.character(equation)){
			stop("cannot coerce equation to a string")}
		}
	if(length(equation)!=1){stop("equation is not a string of length 1")}
	#get arg names
	argNames<-all.vars(as.formula(paste0(paste0(notName,"~"),equation)))
	argNames<-argNames[!(argNames==notName)]
	#create blanks for arguments
	args<-rep(list(bquote()),length(argNames))
	names(args)<-argNames
	#now create an empty function
	newFun<-function(){}
      formals(newFun)<-args
      body(newFun)<-parse(text=equation)
      environment(newFun)<-envir
	return(newFun)
	}
