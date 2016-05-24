#' Joint Probability of A Clade Surviving Infinitely or Being Sampled Once
#'
#' Given the rates of branching, extinction and sampling, calculates the joint
#' probability of a random clade (of unknown size, from 1 to infinite) either
#' (a) never going extinct on an infinite time-scale or (b) being sampled at
#' least once, if it does ever go extinct. As we often assume perfect or close
#' to perfect sampling at the modern (and thus we can blanket assume that living
#' groups are sampled), we refer to this value as the Probability of Being Sampled,
#' or simply P(s). This quantity is useful for calculating the
#' probability distributions of waiting times that depend on a clade being
#' sampled (or not).

#' @details
#' Note that the use of the word 'clade' here can mean a monophyletic group
#' of any size, including a single 'species' (i.e. a single phylogenetic branch)
#' that goes extinct before producing any descendants. Many scientists I have
#' met reserve the word 'clade' for only groups that contain at least one
#' branching event, and thus contain two 'species'. I personally prefer to
#' use the generic term 'lineage' to refer to monophyletic groups of one to
#' infinity members, but others reserve this term for a set of morphospecies
#' that reflect an unbroken anagenetic chain.
#'
#' Obviously the equation used makes assumptions about prior knowledge of the
#' time-scales associated with clades being extant or not: if we're talking
#' about clades that originated a short time before the recent, the clades that
#' will go extinct on an infinite time-scale probably haven't had enough time
#' to actually go extinct. On reasonably long time-scales, however, this infinite
#' assumption should be reasonable approximation, as clades that survive 'forever'
#' in a homogenous birth-death scenario are those that get very large immediately
#' (similarly, most clades that go extinct also go extinct very shortly after
#' originating... yes, life is tough).
#'
#' Both an exact and inexact (iterative) solution is offered; the exact solution
#' was derived in an entirely different fashion but seems to faithfully reproduce
#' the results of the inexact solution and is much faster. Thus, the exact
#' solution is the default. As it would be very simple for any user to look this up
#' in the code anyway, here's the unpublished equation for the exact solution:
#'
#' \eqn{Ps = 1-(((p+q+r)-(sqrt(((p+q+r)^2)-(4*p*q))))/(2*p))}

#' @inheritParams SamplingConv

#' @param useExact If TRUE, an exact solution developed by Emily King is
#' used; if FALSE, an iterative, inexact solution is used, which is somewhat slower
#' (in addition to being inexact...).

#' @return
#' Returns a single numerical value, representing the joint probability of a clade
#' generated under these rates either never going extinct or being sampled before
#' it goes extinct.

#' @seealso
#' \code{\link{SamplingConv}}

#' @author
#' This function is entirely the product of a joint effort between the package author
#' (David W. Bapst), Emily A. King and Matthew W. Pennell. In particular, Emily King
#' solved a nasty bit of calculus to get the inexact solution and later re-derived
#' the function with a quadratic methodology to get the exact solution. Some elements
#' of the underlying random walk model were provided by S. Nalayanan (a user on the
#' website stackexchange.com) who assisted with a handy bit of math involving Catalan
#' numbers.

#' @references
#' Bapst, D. W., E. A. King and M. W. Pennell. In prep. Probability models
#' for branch lengths of paleontological phylogenies.
#' 
#' Bapst, D. W. 2013. A stochastic rate-calibrated method for time-scaling
#' phylogenies of fossil taxa. \emph{Methods in Ecology and Evolution}.
#' 4(8):724-733.
#'

#' @examples
#' #with exact solution
#' pqr2Ps(p=0.1,q=0.1,r=0.1,useExact=TRUE)
#' #with inexact solution
#' pqr2Ps(p=0.1,q=0.1,r=0.1,useExact=TRUE)

#' @export
pqr2Ps<-function(p,q,r,useExact=TRUE){
	#the probability of sampling at least once an extinct clade of unknown size
	if(useExact){
		#test emily's alternative derivation 10-03-13
		#quadratic solution with minus
		res<-((p+q+r)^2)-(4*p*q)
		res<-1-(((p+q+r)-(sqrt(res)))/(2*p))
	}else{
		#original, needs to be iterated over N
		res<-numeric()
		for(N in 1:1000){
			res1<-((p^(N-1))*(q^N)*choose((2*N)-2,N-1))/(N*((p+q+r)^((2*N)-1)))
			if(is.na(res1)){break}	
			if(res1==Inf){break}
			if(res1==0){break}
			res[N]<-res1
			}
		res<-1-sum(res)
		}
	names(res)<-NULL
	return(res)
	}