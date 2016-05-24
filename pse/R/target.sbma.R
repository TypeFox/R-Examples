#' Adaptative generation of Latin Hypercubes
#'
#'	Generates a series of Latin Hypercube Samples for a model until a pair of LHS
#'	present a measure of agreement equal to or greater than a specified target.
#' @inheritParams LHS
#' @param target The desired SBMA.
#' @param init  The size of the initial LHS generated.
#' @param inc The increment between successive runs. For example, if init = 5 and inc = 20, the first
#'	  LHS will be generated with size 5, the second with size 25.
#' @param FUN
#'	  When the model returns more than one response, SBMA values are calculated for each variable.
#'	  The FUN argument specifies how to combine these SBMA values. The recommended default is to
#'	  chose the minimum value.
#' @return Returns the largest LHS generated.
#' @export
target.sbma <- function(target, model, factors,  q = NULL, q.arg = NULL, res.names=NULL, method=c("HL", "random"), 
						opts=list(), init=length(factors)+2, inc=100, FUN=min) {
	#initial LHS
	N = init
	method=match.arg(method)
	print("INFO: initial run...")
	oldL <- LHS(model, factors, N, q, q.arg, res.names, method, opts, nboot=0)
	while (TRUE) {
		N = N + inc
		print(paste("INFO: LHS with N =", N));
		newL <- LHS(model, factors, N, q, q.arg, res.names, method, opts, nboot=0)
		s <- FUN(sbma(newL, oldL))
		print(paste("sbma of ", round(s,3)," (target ",target,")", sep=""))
		if (s >= target) return (newL);
		oldL <- newL;
	}
}

