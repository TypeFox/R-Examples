#' Utility extract functions
#' 
#' Several functions have been added to extract components from externally saved crm models
#' function.wrapper accepts a character string for a model, loads it and then runs the specified
#' function on it. Currently, only 2 functions have been defined: fx.aic to compute aic and
#' fx.par.count to extract parameter count. Currently parameters other than x (eg chat) are passed through
#' the environment. Possibly could have used ...
#' 
#' @aliases function.wrapper fx.aic fx.par.count
#' @usage function.wrapper(x,fx,base="",...)
#' 
#'        fx.aic(x) 
#' 
#'        fx.par.count(x)
#' 
#' @param x character string for the model stored as external object (.rda)
#' @param fx function to be called for each model
#' @param base base name for models
#' @param ... additional values that are added to environment of fx (eg chat for fx.aic)
#' @export function.wrapper
#' @export fx.aic
#' @export fx.par.count
#' @return extracted value defined by function
#' @author Jeff Laake
#' @keywords utility
function.wrapper <- function(x,fx,base="",...)
{
	model.name=paste(x,collapse=".")
	e1=new.env()
	dots=list(...)
	sapply(1:length(dots),function(x) assign(names(dots)[x],dots[[x]],envir=e1))
	environment(fx)=e1
	eval(parse(text=paste('load(file="',base,model.name,'.rda")',sep="")),envir=e1)
	eval(parse(text=paste("result<-fx(",model.name,")",sep="")),envir=e1)	
	return(get("result",envir=e1))
}

fx.aic=function(x) 
{
	if(!"chat"%in% ls())chat=1
	x$neg2lnl/chat + 2*length(x$beta)
}
fx.par.count=function(x) length(grep(paste("^",par,":",sep=""),names(x$beta)))
