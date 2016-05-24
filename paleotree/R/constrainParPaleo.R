#' Constrain Parameters for a Model Function from paleotree
#'
#' This function constrains a model to make submodels with fewer parameters,
#' using a structure and syntax taken from the function \code{constrain}
#' in Rich Fitzjohn's package \code{diversitree}.

#' @details
#' This function is based on (but does not depend on) the function \code{constrain}
#' from the package \code{diversitree}. Users should refer to this parent function for
#' more detailed discussion of model constraint usage and detailed examples.
#'
#' The parent function was forked to add functionality necessary for dealing with the
#' high parameter count models typical to some paleontological analyses, particularly
#' the inverse survivorship method. This necessitated that the new function be entirely
#' separate from its parent. Names of functions involved (both exported and not)
#' have been altered to avoid overlap in the package namespaces. Going forward,
#' the paleotree package maintainer (Bapst) will try to be vigilant with
#' respect to changes in \code{constrain} in the original package, \code{diversitree}.
#'
#' Useful information from the diversitree manual (11/01/13):
#'
#' "If \code{f} is a function that takes a vector \code{x} as its first
#' argument, this function returns a new function that takes a
#' shorter vector \code{x} with some elements constrained in some way;
#' parameters can be fixed to particular values, constrained to be the
#' same as other parameters, or arbitrary expressions of free
#' parameters."
#'
#' In general, formulae should be of the structure:
#' 
#' \emph{LHS ~ RHS}
#'
#' ...where the LHS is the 'Parameter We Want to Constrain' and the
#' RHS is whatever we are constraining the LHS to, usually another
#' parameter. LHS and RHS are the 'left-hand side' and
#' 'right-hand side' respectively (which I personally find obscure).
#'
#' Like the original \code{constrain} function this function is based on,
#' this function cannot remove constraints previously placed on a model
#' object and there may be cases in which the constrained function may not 
#' make sense, leading to an error. The original function will sometimes
#' issue nonsensical functions with an incorrect number/names of parameters
#' if the parameters to be constrained are given in the wrong order in
#' formulae.
#'
#' \subsection{Differences from diversitree's constrain Function}{
#'
#' This forked paleotree version of constrain has two additional features,
#' both introduced to aid in constraining models with a high number of 
#' repetitive parameters. (I didn't invent these models, don't shoot the messenger.)
#'
#' First, it allows nuanced control over the constraining of many
#' parameters simultaneously, using the 'all' and 'match' descriptors. This
#' system depends on parameters being named as such: name.group1.group2.group3
#' and so on. Each 'group' is reference to a system of groups, perhaps referring to a
#' time interval, region, morphology, taxonomic group or some other discrete
#' characterization among the data (almost all functions envisioned for
#' paleotree are for per-taxon diversification models). The number of group systems
#' is arbitrary, and may be from zero to a very large number; it depends on the
#' 'make' function used and the arguments selected by the user. For example, the
#' parameter 'x.1' would be for the parameter 'x' in the first group of the first group
#' system (generally a time interval for most paleotree functions). For a more
#' complicated exampled, with the parameter 'x.1.3.1', the third group
#' for the second group system (perhaps this taxonomic data point has a morphological
#' feature not seen in some other taxa) and group 1 of the third group system (maybe
#' biogeographic region 1? the possibilities are endless depending on user choices). 
#'
#' The 'all' option work like so: if 'x.all~x.1' is given as a formulae, then all x
#' parameters will be constrained to equal x.1. For example, if there is x.1, x.2,
#' x.3 and x.4 parameters for a model, 'x.all~x.1' would be equivalent to
#' individually giving the formulae 'x.2~x.1', 'x.3~x.1' and 'x.4~x.1'. This
#' means that if there are many parameters of a particular type (say, 50 'x'
#' parameters) it is easy to constrain all with a short expression. It is not
#' necessary that the  The 'all' can be used anywhere in the name of the parameter
#' in a formulae, including to make all parameters for a given 'group' the same.
#' Furthermore, the LHS and RHS don't need to be same parameter group, and both can
#' contain 'all' statements, even \emph{multiple} 'all' statements. Consider these
#' examples, all of which are legal: 
#'
#' \describe{ 
#'  \item{x.all ~ y.1}{Constrains all values of the parameter x for every group to be
#' equal to the single value for the parameter y for group 1 (note that there's only
#' a single set of groups).}

#'  \item{all.1 ~ x.1}{Constrains all parameters for the first group to equal each other, 
#' here arbitrary named x.1. For example, if there is parameters named x.1, y.1 and z.1,
#' all will be constrained to be equal to a single parameter value.}

#'  \item{x.all.all ~ y.2.3}{Constrains all values for x in any and all groups to
#' equal the single value for y for group 2 (of system 1)  and group 3 (of system 2).}

#'  \item{x.all ~ y.all}{Constrains all values of x for every group and y for every
#' group to be equal to a \emph{single} value, which by default will be reported as y.1}}
#'
#' The 'match' term is similar, allowing parameter values from the same group
#' to be quickly matched and made equivalent. These 'match' terms must have a
#' matching (hah!) term both in the corresponding LHS and RHS of the formula.
#' For example, consider 'x.match~y.match' where there are six parameters: x.1,
#' x.2, x.3, y.1, y.2 and y.3. This will effectively constrain x.1~y.1, x.2~y.2
#' and x.3~y.3. This is efficient for cases where we have some parameters that
#' we often treat as equal. For example, in paleontology, we sometimes make a
#' simplifying assumption that birth and death rates are equal in multiple
#' time intervals. Some additional legal examples are:
#'
#' \describe{ 
#' \item{x.match.1 ~ y.match.1}{This will constrain only parameters of x and y to
#' to equal each other if they both belong to the same group for the first group
#' system AND belong to group 1 of the first group.}

#' \item{all.match. ~ x.match}{This will constrain all named parameters in each
#' group to equal each other; for example, if there are parameters x.1, y.1, z.1,
#' x.2, y.2 and z.2, this will constrain them such that y.1~x.1, z.1~x.1, y.2~x.2
#' and z.2~x.2, leaving x.1 and x.2 as the only parameters effectively.}}
#'
#' There are two less fortunate qualities to the introduction of the above terminology.
#'
#' Unfortunately, this package author apologizes that his programming skills are
#' not good enough to allow more complex sets of constraints, as would be typical
#' with the parent function, when 'all' or 'match' terms are included. For example,
#' it would not be legal to attempt to constraint 'y.all ~ x.1 / 2', where the user
#' presumably is attempting to constrain all y values to equal the x parameter
#' to equal half of the x parameter for group 1. This will not be parsed as such
#' and should return an error. However, there are workarounds, but they require
#' using constrainParPaleo more than once. For the above example, a user could
#' first use 'y.all ~ y.1' constraining all y values to be equal. Then a user
#' could constrain with the formula 'y.1 ~ x.1 / 2' which would then constrain
#' y.1 (and all the y values constrained to equal it) to be equal to the desired
#' fraction.
#'
#' Furthermore, this function expects that parameter names don't already have
#' period-separated terms that are identical to 'all' or 'match'. No function 
#' in paleotree should produce such natively. If such were to occur, perhaps
#' by specially replacing parameter names, constrainParPaleo would confuse
#' these terms for the specialty terms described here.
#' 
#' Secondly, this altered version of constrain handles the parameter bounds included as
#' attributes in functions output by the various 'make' functions. This means that if
#' 'x.1 ~ y.1', constrainParPaleo will test if the bounds on x.1 and y.1 are the same.
#' If the bounds are not the same, constrainParPaleo will return an error.
#' This is important, as some models in paleotree may make a parameter a rate (bounded
#' zero to some value greater than one) or a probability (bounded zero to one),
#' depending on user arguments. Users may not realize these differences and, in many
#' cases, constraining a rate to equal a probability is nonsense (absolute poppycock!).
#' If a user really wishes to constrain two parameters with different bounds to be equal
#' (I have no idea why anyone would want to do this), they can use the parameter bound
#' replacement functions described in \code{\link{modelMethods}} to set the parameter
#' bounds as equal. Finally, once parameters with the same bounds are constrained, the
#' output has updated bounds that reflect the new set of parameters
#' for the new constrained function.}
#'
  
#' @param f A function to constrain. This function must be of S3 class
#' 'paleotreeFunc' and have all necessary attributes expected of that
#' class, which include parameter names and upper and lower bounds. As
#' I have deliberately not exported the function which creates this class,
#' it should be impossible for regular users to obtain such objects easily
#' without using one of the 'make' functions, which automatically output
#' a function of the appropriate class and attributes.

#' @param ... Formulae indicating how the function should be constrained.
#' See details and examples for lengthy discussion.

#' @param formulae Optional list of constraints, possibly in addition to
#' those in \code{...}

#' @param names Optional Character vector of names, the same length as
#' the number of parameters in \code{x}.  Use this only if
#' \code{\link{parnames}} does not return a vector for your function.
#' Generally this should not be used. DWB: This argument is kept for
#' purposes of keeping the function as close to the parent as possible
#' but, in general, should not be used because the input function must
#' have all attributes expected of class 'paleotreeFunc', including
#' parameter names.

# @param bounds Optional list composed of two numerical vectors of the
# same length as the number of parameters, representing the bounds on
# the parameters, if those are not given by your function; i.e. 
# \code{\link{parbounds}} does not return an apppropriate list for your
# function.

#' @param extra Optional vector of additional names that might appear on
#' the RHS of constraints but do not represent names in the function's
#' \code{argnames}.  This can be used to set up dummy variables
#' (example coming later).

#' @return 
#' Modified from the diversitree manual:

#' This function returns a constrained function that can be passed
#'  through to the optimization functions of a user's choice, such as
#' \code{\link{optim}}, \code{find.mle} in diversitree or \code{mcmc}.
#' It will behave like any other function.  However, it has a modified
#' \code{class} attribute so that some methods will dispatch differently:
#' \code{\link{parnames}}, for example, will return the names of the
#' parameters of the constrained function and \code{\link{parInit}} will
#' return the initial values for those same constrained set of parameters.
#' All arguments in addition to \code{x} will be passed through to the
#' original function \code{f}.
#'
#' Additional useful information from the diversitree manual (11/01/13):
#'
#' For help in designing constrained models, the returned function has
#' an additional argument \code{pars.only}, when this is \code{TRUE} the
#' function will return a named vector of arguments rather than evaluate
#' the function (see Examples).

#' @seealso
#' As noted above, this function is based on (but does not depend on) the
#' function \code{constrain} from the library \code{diversitree}.

#' @author 
#' This function (and even this help file!) was originally written by Rich
#' Fitzjohn for his library \code{diversitree}, and subsequently rewritten
#' and modified by David Bapst.

#' @references
#' FitzJohn, R. G. 2012. Diversitree: comparative phylogenetic analyses of
#' diversification in R. \emph{Methods in Ecology and Evolution} 3(6):1084-1092.

#' @examples
#' #simulation example with make_durationFreqCont, with three random groups
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' grp1 <- matrix(sample(1:3,nrow(taxa),replace=TRUE),,1)   #groupings matrix
#' likFun <- make_durationFreqCont(rangesCont,groups=grp1)
#' 
#' # can constrain both extinction rates to be equal
#' constrainFun <- constrainParPaleo(likFun,q.2~q.1)
#' 
#' #see the change in parameter names and bounds
#' parnames(likFun)
#' parnames(constrainFun)
#' parbounds(likFun)
#' parbounds(constrainFun)
#' 
#' # some more ways to constrain stuff!
#' 
#' #constrain all extinction rates to be equal
#' constrainFun <- constrainParPaleo(likFun,q.all~q.1)
#' parnames(constrainFun)
#' 
#' #constrain all rates for everything to be a single parameter
#' constrainFun <- constrainParPaleo(likFun,r.all~q.all)
#' parnames(constrainFun)
#' 
#' #constrain all extinction rates to be equal & all sampling to be equal
#' constrainFun <- constrainParPaleo(likFun,q.all~q.1,r.all~r.1)
#' parnames(constrainFun)
#' 
#' #similarly, can use match.all to make all matching parameters equal each other
#' constrainFun <- constrainParPaleo(likFun,match.all~match.all)
#' parnames(constrainFun)
#' 
#' #Constrain rates in same group to be equal
#' constrainFun <- constrainParPaleo(likFun,r.match~q.match)
#' parnames(constrainFun)
#' 

#' @export
constrainParPaleo<-function(f, ..., formulae=NULL, names=parnames(f),extra=NULL) {
	#based on Rich FitzJohn's constrain function for diversitree 10-22-13
		# comment lines with double ## indicate Rich's original comments
		#I claim all and any responsibility for how fugly I can make Rich's code
			#...its actually kind of a challenge - DWB
	##
	## For the first case, everything is OK on the lhs and rhs
	## For subsequent cases:
	## lhs cannot contain things that are
	## - constrained things (already lhs anywhere)
	## - things constrained to (things on the rhs anywhere)
	## rhs cannot contain things that are
	## - constrained things (already lhs anywhere)
	## It is possibly worth pulling out all the numerical constants and
	## the "paired" parameters here to avoid using eval where
	## unnecessary. However, this makes the function substantially uglier
	## for a very minor speedup.
	#
	#let's make some example data
		#f=function(pqr=c(p,q,r)){pqr[1]^2+2*pqr[2]+pqr[3]} 
		#f<-make_paleotreeFunc(f,c("p","q","r"),list(c(0,0,0),rep(Inf,3)))
		#formulae=c(p~q,list());names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
		#f=function(pqr=c(p,q,r)){p^2+2*q+r}; 
		#f<-make_paleotreeFunc(f,c("p","q","r"),list(c(0,0,0),rep(Inf,3)))
		#formulae=c(p~q,r~q,list());names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
		#f=function(pqr=c(p.1,p.2,q.1,q.2,r.1,r.2)){p.1^2+2*q.1+r.1/(p.2^2+2*q.2+r.2)}
		#f<-make_paleotreeFunc(f,c("p.1","q.1","r.1","p.2","q.2","r.2"),list(rep(0,6),rep(Inf,6)))
		#formulae=c(p.1~q.all,p.match~r.match,list());names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
		#FOR USE WITH known f
		#formulae = c(rRate~pRate);names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
	#get bounds
	bounds <- parbounds(f)
	#
	if(!is(f,'paleotreeFunc')){
		stop("Given function does not appear to be a paleotree likelihood function")}
	if ( inherits(f, "constrained") ) {	#this thing checks to see if its already constrained 
		formulae <- c(attr(f, "formulae"), formulae)
		f <- attr(f, "origfunction")
		}
	rels <- list()				#rels are the things we're gonna constrain to something else
	names.lhs <- names.rhs <- names	#lhs is untouched pars??, rhs is the final pars??
	formulae <- c(formulae, list(...))	#adding the ... to formulae
	#expand formulae in case they contain systematic constraints here! here's some examples:
		#names=c("p.1.1","q.1.1","p.2.1","q.2.1","p.1.2","q.1.2","p.2.2","q.2.2")
		#formulae = c(p.all.match~q.all.match,list())
		#formulae = c(p.1.1~p.all.all,list())
		#formulae = c(p.1.match~q.all.match,list())
		#formulae = c(p.1.match~q.all.match,p.1.1~p.all.all,list())
	breakTerms<-lapply(formulae,function(x) unlist(strsplit(all.vars(as.formula(x)),".",fixed=TRUE)))
	needExpand<-sapply(breakTerms,function(x) any("match"==x)|any("all"==x))
	if(any(needExpand)){
		breakNames<-t(rbind(sapply(names,function(x) unlist(strsplit(x,".",fixed=TRUE)))))
		nparcat<-ncol(breakNames)
		newFormulae<-list()
		for(i in which(needExpand)){
			newFormulae<-c(newFormulae,expandConstrainForm(formula=formulae[[i]],breakNames=breakNames,nparcat=nparcat))
			}
		formulae<-c(formulae[-which(needExpand)],newFormulae) #formulae		
		formulae<-unique(formulae)	#any duplicates?
		}
	for( formula in formulae ) {
		res <- constrainParsePaleo(formula, names.lhs, names.rhs, extra)
		if ( attr(res, "lhs.is.target") ) {
			i <- try(which( sapply(rels,function(x) identical(x, res[[1]]))),silent=TRUE)
			if(inherits(i,"try-error")){
				stop(sprintf("Error parsing constraint with %s on lhs",as.character(res[[1]])))
				}
			rels[i] <- res[[2]]	#DWB: gives warning message that symbol cannot be coerced to list
			## This will not work with *expressions* involving the LHS; that
			## would require rewriting the expressions themselves (which
			## would not be too hard to do). But for now let's just cause
			## an error...
			lhs.txt <- as.character(res[[1]])
			if ( any(sapply(rels, function(x) lhs.txt %in% all.vars(x))) ){
				stop(sprintf("lhs (%s) is in an expression and can't be constrained",lhs.txt))
				}
			}
		names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
		names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
		rels <- c(rels, structure(res[2], names=as.character(res[[1]])))
	  	}
	#in a function (p,q,r), with constraint p~q, names.lhs="r", names.rhs="q""r", rels = list(p=r)
	#okay, now we know which ones will be lhs, rhs and rels
	#need to test that all the bounds for the equivalencies are the same
		#check each rels for consistency with bounds
	relsIsPar<-rels[sapply(rels,function(x) names==x)]	#test to see which are pars
	if(length(relsIsPar)>0){
		termBound<-t(sapply(names(relsIsPar),function(x) 
			sapply(bounds,function(y) y[which(names==x)])))
		relsBound<-t(sapply(relsIsPar,function(x) 
			sapply(bounds,function(y) y[which(names==x)])))
		colnames(relsBound)<-colnames(termBound)<-NULL
		if(!identical(relsBound,termBound)){
			noMatch<-which(!apply(termBound==relsBound,1,all))
			noMatch<-paste(i,"~",names(relsIsPar)[noMatch])
			stop(paste("Upper and Lower bounds do not match for",noMatch))
			}
		}
	#back to usual diversitree code for a moment
	i <- match(unique(sapply(rels, as.character)), extra)	#match
	final <- c(extra[sort(i[!is.na(i)])], names.rhs)
	npar <- length(final)	
	#need to update the bounds at the same time the pars get updated
	bounds<-lapply(bounds,function(x) x[sapply(final,function(x) which(x==names))])
	## "free" are the parameters that have nothing special on their RHS
	## and are therefore passed directly through
	free <- setdiff(names.rhs, names(rels))
	free.i <- match(free, names) # index in full variables
	free.j <- match(free, final) # index in given variables.
	## Targets are processed in the same order as given by formulae.
	target.i <- match(names(rels), names)
	pars.out <- rep(NA, length(names))
	names(pars.out) <- names
	g <- function(pars, ..., pars.only=FALSE) {
		if ( length(pars) != npar ){
			stop(sprintf("Incorrect parameter length: expected %d, got %d",npar, length(pars)))
			}
	    pars.out[free.i] <- pars[free.j]
	    e <- structure(as.list(pars), names=final)
	    pars.out[target.i] <- unlist(lapply(rels, eval, e))
		if(pars.only){
			res<-pars.out
		}else{
			res<-f(pars.out, ...)
			}
		return(res)
		}
	class(g) <- c("constrained", class(f))
	attr(g, "parnames") <- final
	attr(g, "parbounds") <- bounds
	attr(g, "formulae") <- formulae
	attr(g, "extra") <- extra
	attr(g, "origfunction") <- f
	return(g)
	}