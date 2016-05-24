#' Model Function Methods: Parameter Names, Bounds and Initial Values
#'
#' A large number of functions for obtaining and modifying the parameters
#' of likelihood models made in paleotree. These functions allow users to obtain
#' or set parameter names, or obtain and set parameter bounds, both of which
#' are treated as an attribute of the function class used by paleotree. In
#' practice, this allows users to quickly obtain parameter names and upper
#' and lower values for use in bounded optimizers, including reasonable
#' starting values.

#' @details
#' Parameter names cannot be changed for a constrained function.
#'
#' The parInit function calls the bounds for each parameter and gives a randomly
#' selected value selected from a uniform distribution, using the parameter bounds
#' for each parameter as the bounds on the uniform distribution. This users a
#' shorthand to quickly generate initial parameter values which are within the
#' set bounds, for use in functions such as \code{\link{optim}}. The random
#' sampling of initial values allows a user to quickly assess if initial
#' parameter values affect the optimization by simply rerunning the function on new values.
#' Infinite initial parameter values (resulting from infinite bounds) are discarded, and
#' replaced with the lower bound value (assuming only upper bounds are infinite...).
#' Some randomly selected initial parameter values may be too high (due to the liberal
#' upper bounds I set for parameters in many of the likelihood functions) and
#' thus users should always try slightly different values to see if the resulting
#' maximum likelihood parameter values change.
#'
#' As parInit depends on the upper and lower bounds attribute, no function is offered
#' to allow it to be replaced (as there is nothing to replace!).

#' @param x A function of S3 class 'paleotreeFunc' with all necessary attributes
#' expected of that class, which include parameter names and upper and lower bounds.
#' As I have deliberately not exported the function which creates this class, it
#' should be impossible for regular users to obtain such objects easily without
#' using one of the 'make' functions, which automatically output a function of the
#' appropriate class and attributes.

#' @param ... 'Ignored arguments to future methods' (i.e. for diversitree). Kept here only
#' so constrainParPaleo is kept as close to the parent method in diversitree as possible.

#' @param value The new value with which to replace the parameter names or bounds. Must
#' be a vector of the same length as the number of parameters. For parbounds, must
#' be a list composed of two vectors.

#' @return
#' Returns the sought parameter names, bounds or initial values or (for the replacement methods) 
#' returns a modified function with the respective attributes altered.

#' @aliases modelMethods parbounds parbounds.constrained parbounds.paleotreeFunc
#' parInit parInit.constrained parInit.paleotreeFunc parLower
#' parLower.constrained parLower.paleotreeFunc parnames
#' parnames.constrained parnames.paleotreeFunc parUpper
#' parUpper.constrained parUpper.paleotreeFunc

#' @seealso
#' These model methods were introduced to interact with the new model framework introduced in
#' paleotree v1.9, in particular to interface with \code{\link{constrainParPaleo}}.

#' @author
#' These functions are strongly based on or inspired by the \code{argnames} functions
#' provided for handling models in Rich Fitzjohn's library \code{diversitree}, but
#' the functions presented here are derviations written by David Bapst.

#' @examples
#' #example with make_durationFreqCont
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' likFun <- make_durationFreqCont(rangesCont)
#' 
#' #get parameter names
#' parnames(likFun)
#' 
#' #get the bounds for those parameters
#' parbounds(likFun)
#' 
#' #can also get these seperately
#' parLower(likFun)
#' parUpper(likFun)
#' 
#' #initial parameter values
#' parInit(likFun)   #arbitrary midway value between par bounds
#' 
#' #can then use these in optimizers, such as optim with L-BFGS-B
#' #see the example for make_durationFreqCont
#' 
#' #renaming parameter names
#' likFun2 <- likFun
#' parnames(likFun2) <- c("extRate","sampRate")
#' parnames(likFun2)
#' #test if reset correctly
#' parnames(likFun2)==c("extRate","sampRate")

#' #also works for constrained functions
#' constrainFun<-constrainParPaleo(likFun,q.1~r.1)
#' parnames(constrainFun)
#' #also modified the parameter bounds, see!
#' parbounds(constrainFun)
#' parInit(constrainFun)
#' #but cannot rename parameter for constrained function!
#'

#parnames

#' @name modelMethods
#' @rdname modelMethods
#' @export
parnames <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parnames")
	}

#' @rdname modelMethods
#' @export
parnames.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parnames")
	}

#' @rdname modelMethods
#' @export
parnames.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parnames")
	}
	
#' @rdname modelMethods
#' @export
`parnames<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parnames<-")
	}
	
#' @rdname modelMethods
#' @export
`parnames<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parnames on a constrained function")
	}
	
#' @rdname modelMethods
#' @export
`parnames<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	#parnames<-attr(x,"parnames")
	value<-as.character(value)
	if(length(value)!=np){stop("length of new parnames not equal to number of parameters")}
	if(any(is.na(value))){stop("NA values in parnames replacement")}
	if(any(duplicated(value))){stop("Duplicated names in parnames replacement")}
	attr(x,"parnames")<-value
	return(x)
	}

#parbounds

#' @rdname modelMethods
#' @export
parbounds <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parbounds")
	}

#' @rdname modelMethods	
#' @export
parbounds.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parbounds")
	}
	
#' @rdname modelMethods
#' @export
parbounds.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parbounds")
	}
	
#' @rdname modelMethods
#' @export
`parbounds<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parbound<-")
	}
	
#' @rdname modelMethods
#' @export
`parbounds<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parbounds on a constrained function")
	}
	
#' @rdname modelMethods
#' @export
`parbounds<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	if(!is.list(value) | !length(value)==2){stop("parbounds needs to be a list composed of two vectors")}
	lower<-as.numeric(value[[1]])
	if(length(lower)!=np){stop("length of new lower parbounds not equal to number of parameters")}
	if(any(is.na(lower))){stop("NA values in lower parbounds replacement")}
	upper<-as.numeric(value[[2]])
	if(length(upper)!=np){stop("length of new upper parbounds not equal to number of parameters")}
	if(any(is.na(upper))){stop("NA values in upper parbounds replacement")}
	attr(x,"parbounds")<-value
	return(x)
	}
	

#get upper and lower bounds for parameters

#' @rdname modelMethods
#' @export
parLower <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parLower")
	}
	
#' @rdname modelMethods
#' @export
parLower.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parbounds")[[1]]
	}
	
#' @rdname modelMethods
#' @export
parLower.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parbounds")[[1]]
	}
	
#' @rdname modelMethods
#' @export
`parLower<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parLower<-")
	}
	
#' @rdname modelMethods
#' @export
`parLower<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parLower on a constrained function")
	}
	
#' @rdname modelMethods
#' @export
`parLower<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	lower<-as.numeric(value)
	if(length(lower)!=np){stop("length of new lower parbounds not equal to number of parameters")}
	if(any(is.na(lower))){stop("NA values in lower parbounds replacement")}
	attr(x,"parbounds")[[1]]<-value
	return(x)
	}
	
#' @rdname modelMethods
#' @export
parUpper <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parUpper")
	}
	
#' @rdname modelMethods
#' @export
parUpper.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parbounds")[[2]]
	}
	
#' @rdname modelMethods
#' @export
parUpper.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parbounds")[[2]]
	}
	
#' @rdname modelMethods
#' @export
`parUpper<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parUpper<-")
	}
	
#' @rdname modelMethods
#' @export
`parUpper<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parUpper on a constrained function")
	}
	
#' @rdname modelMethods
#' @export
`parUpper<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	upper<-as.numeric(value)
	if(length(upper)!=np){stop("length of new upper parbounds not equal to number of parameters")}
	if(any(is.na(upper))){stop("NA values in upper parbounds replacement")}
	attr(x,"parbounds")[[2]]<-value
	return(x)
	}

#initial parameter values

#' @rdname modelMethods
#' @export
parInit <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parInit")
	}

#' @rdname modelMethods
#' @export
parInit.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	#res<-(attr(x, "parbounds")[[2]]-attr(x, "parbounds")[[1]])/2
	res<-runif(length(attr(x,"parnames")),attr(x, "parbounds")[[1]],attr(x, "parbounds")[[2]])
	#infinite bounds probably too far from actual param value; use lower bound instead
	if(any(is.infinite(res))){res<-attr(x, "parbounds")[[1]]}
	names(res)<-attr(x,"parnames")
	return(res)
	}

#' @rdname modelMethods	
#' @export
parInit.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	#res<-(attr(x, "parbounds")[[2]]-attr(x, "parbounds")[[1]])/2
	res<-runif(length(attr(x,"parnames")),attr(x, "parbounds")[[1]],attr(x, "parbounds")[[2]])
	#infinite bounds probably too far from actual param value; use lower bound instead
	if(any(is.infinite(res))){res<-attr(x, "parbounds")[[1]]}
	names(res)<-attr(x,"parnames")
	return(res)
	}