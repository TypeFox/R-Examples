#' Dual optimization for Kiefer-Wolfowitz problems
#' 
#' Interface function for calls to optimizer from various REBayes functions
#' There are currently two options for the optimization:  Mosek (the default)
#' is the original, preferred option and uses interior point methods.
#' It relies on the \pkg{Rmosek} interface to R available from 
#' \url{http://rmosek.r-forge.r-project.org} 
#' See userguide for Rmosek for further details.  A more experimental option
#' employs the \pkg{pogs} package available from \url{https://github.com/foges/pogs}
#' and employs an ADMM (Alternating Direction Method of Multipliers) approach.
#' 
#' @param A Linear constraint matrix
#' @param d constraint vector
#' @param w weights for \code{x} should sum to one.
#' @param ...  other parameters passed to control optimization:  These may
#' include \code{rtol} the relative tolerance for dual gap convergence criterion,
#' \code{verb} to control verbosity desired from mosek, \code{verb = 0} is quiet,
#' \code{verb = 5} produces a fairly detailed iteration log,
#' \code{method} controls the choice of optimizer:  by default this is "mosek"
#' which employs interior point methods, however if \code{method = "pogs"} then
#' optimization is carried out by the ADMM methods described in Fougner and
#' Boyd (2015).  This is a first order descent method most suitable for large
#' problems for which parallelization is desirable.  For most REBayes applications
#' the default "mosek" method is appropriate and "pogs" should be considered
#' experimental.  Note that there is not yet a "pogs" implementation for \code{medde} problems.
#' Note also that \code{method = "pogs"} assumes a distinct control list.
#' Users are responsible for specifying correctly named control variables for each method.
#' The most advantageous implementation of "pogs" requires (CUDA) GPU hardware.
#' \code{control} is a control list consisting of sublists \code{iparam},
#' \code{dparam}, and \code{sparam}, containing elements of various mosek
#' control parameters.  See the Rmosek and Mosek manuals for further details.
#' A prime example is \code{rtol} which should eventually be deprecated and
#' folded into \code{control}, but will persist for a while for compatibility
#' reasons.  The default for \code{rtol} is 1e-6, but in some cases it is
#' desirable to tighten this, say to 1e-10.  Another example that motivated the introduction of
#' \code{control} would be \code{control = list(iparam = list(num_threads =
#' 1))}, which forces Mosek to use a single threaded process.  The default
#' allows Mosek to uses multiple threads (cores) if available, which is
#' generally desirable, but may have unintended (undesirable) consequences when running
#' simulations on clusters.
#' @return Returns a list with components: \item{f}{dual solution vector, the
#' mixing density} \item{g}{primal solution vector, the mixture density
#' evaluated at the data points} \item{logLik}{log likelihood}
#' \item{status}{return status from Mosek}
#' @author R. Koenker
#' @references
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' \emph{JASA}, 109, 674--685.
#'
#' Mosek Aps (2015) Users Guide to the R-to-Mosek Optimization Interface, 
#' \url{http://rmosek.r-forge.r-project.org}.
#' 
#' Fougner, C. (2015) POGS: Proximal Operator Graph Solver, R Package available from
#' \url{http://foges.github.io/pogs}.
#' 
#' Fougner, C. and S. Boyd, (2015) Parameter Selection and Pre-Conditioning for a
#' Graph Form Solver, Stanford Technical Report.
#' @keywords nonparametrics
#' @importFrom methods as new
#' @export
KWDual <- function(A, d, w, ...){
# Dual Kiefer-Wolfowitz MLE for Mixture Problems
#
#       This version implements a class of density estimators solving:
#
#       min_x  {F(x) := sum -log (x_i)}  s.t. A' x <= d, 0 <= x,  
#
#
#	where e.g.  A = phi(outer(Y,g,"fun")), with Y data and g a grid on the support of Y,
#	and "fun"  is some function representing the dependence of the base distribution.
#
#-------------------------------------------------------------------------------------
#
# Roger Koenker 
#
# First version:24 Feb 2012  
# Revised:	10 Jun 2015 # Simplified signature
# Revised:	 2 Jul 2015 # Added pogs method

KWpogs <- function (A, d, w, control) { # POGS implementation of KWDual
    n <- nrow(A)
    m <- ncol(A)
    # Uncomment the next two lines if you want to use pogs
    #f <- list(h = pogs::kIndBox01(n), c = d)
    #g <- list(h = pogs::kNegLog(m), c = w)
    pogs.control <- function(rel_tol=1e-4, abs_tol=1e-4, rho=1.0,
       max_iter=1000, verbose = 1, adaptive_rho=TRUE)
       list(rel_tol=rel_tol, abs_tol=abs_tol, rho=rho,
          max_iter=max_iter, verbose=verbose, adaptive_rho=adaptive_rho)
    params <- pogs.control()
    if(length(control)){
	control <- as.list(control)
	params[names(control)] <- control
    }
    # Uncomment the next line if you want to use pogs
    #z <- pogs::pogs(A, f, g, params)
    # This abuse of notation is needed to conform to KWDual
    f <- z$v/d
    g <- as.vector(t(A) %*% (f * d))
    list(f = f, g = g, status = z$status)
}

n <- nrow(A)
m <- ncol(A)
A <- t(A) 

dots <- list(...)

if(length(dots$method))
    if(dots$method == "pogs") {
	dots$method <- NULL
	return(KWpogs(A, d, w, control = dots))
    }
    else if(!dots$method == "mosek") 
	stop(paste("No applicable KWDual method: ", dots$method))

# Default mosek method
rtol <- ifelse(length(dots$rtol), dots$rtol, 1e-6)
verb <- ifelse(length(dots$verb), dots$verb, 0)
if(length(dots$control)) control <- dots$control
else control <- NULL


C <- rep(0,n)
P <- list(sense = "min")
P$c <- C
P$A <- Matrix::Matrix(A, sparse = TRUE)
P$bc <- rbind(rep(0,m),d)
P$bx <- rbind(rep(0,n),rep(Inf,n))
opro <- matrix ( list (), nrow =5, ncol = n)
rownames ( opro ) <- c(" type ","j","f","g","h")

opro[1,] <-  as.list(rep('log',n))
opro[2,] <-  as.list(1:n)
opro[3,] <-  as.list(-w)
opro[4,] <-  as.list(rep(1,n))
opro[5,] <-  as.list(rep(0,n))
P$scopt<- list(opro = opro)
P$dparam$intpnt_nl_tol_rel_gap <- rtol
if(length(control)){
    P$iparam <- control$iparam
    P$dparam <- control$dparam
    P$sparam <- control$sparam
}
z <- Rmosek::mosek(P, opts = list(verbose = verb))
status <- z$sol$itr$solsta
f <- z$sol$itr$suc
g <- as.vector(t(A) %*% (f * d))
list(f = f, g = g, status = status)
}
