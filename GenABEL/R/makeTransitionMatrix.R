#' Genotype transition probabilities matrices
#' 
#' Function to generate genotypic transition probablilites 
#' matrices, which represent conditional probabilities 
#' P(g1|g2,nmeioses), where g1 is henotype of person one
#' (AA, AB or BB), g2 is genotype of person two, and 
#' nmeioses is the number of meioses separating these 
#' two individuals (0 for twins, 1 for parent-offspring, 
#' c(2,2) for sibs, 2 for grandparent-grandchild pairs, etc.)
#' 
#' @param nmeioses number of meioses separating two 
#' individuals ((a vector) of non-negative integers). 
#' If a vector, it is assumed it lists all meiotic paths
#' connecting the pair
#' @param q (a vector of) the coded allele frequency(ies) 
#' (e.g. "Q.2" of GenABEL-package)
#' 
#' @return If q is scalar, a 3x3 matrix is returned,
#' where elements represent conditional transition 
#' probabilities P(g1|g2,nmeioses); rows correspond to 
#' the genotypes of g1, and columns correspond to the 
#' genotypes of g2. If coded allele is 'B', then 
#' e.g. element [1,2] gives the probability 
#' P(g1='AA'|g2='AB',nmeioses=nmeioses).
#' 
#' If q is a vector, series of above-described matrices 
#' are returned as an 'array' object. A matrix constructed 
#' for certain element q[i] can be accessed via 
#' result[,,i].
#' 
#' @examples
#' # transition matrix for parent-offspring, for q=0.1
#' makeTransitionMatrix(0.1,nmeioses=1)
#' # for a set of q's 
#' makeTransitionMatrix(c(0.1,0.9),nmeioses=1)
#' # for sibs
#' makeTransitionMatrix(0.1,nmeioses=c(2,2))
#' # for half-sibs (or grandparent-grandchild)
#' makeTransitionMatrix(0.1,nmeioses=2)
#' # for remote relatives
#' makeTransitionMatrix(0.1,nmeioses=10)
#' # for independent
#' makeTransitionMatrix(0.1,nmeioses=1000)
#' 
#' @author Yurii Aulchenko
#'   
#' 
makeTransitionMatrix <- function(q,nmeioses=1000) {
# do sanity checks
	if (!is(q,"numeric")) 
		warning("q should be real (scalar or vector) between 0 and 1")
	if (any(q<0) | any(q>1)) 
		stop("(some) q<0 or q>1")
	if (!is(nmeioses,"numeric") || min(nmeioses)<0) 
		stop("nmeioses should be a positive integer or 0")
	if (length(nmeioses)>2) 
		stop("current implementation assumes 2 meiotic paths max")
	# compute transition matrix: P(g1|g2) where g1 is first (rows) and g2 is second dimention (columns)
	getTT <- function(q,pIBD) {
		p <- 1-q
		p2 <- p*p; pq <- p*q; pq2 <- 2*pq; q2 <- q*q; 
		if (length(pIBD)==2) {
			J0 <- (1-pIBD[1])*(1-pIBD[2]); 
			J1 <- pIBD[1]*(1-pIBD[2]) + (1-pIBD[1])*pIBD[2]; 
			J2 <- pIBD[1]*pIBD[2]
		} else {
			J0 <- (1-pIBD) 
			J1 <- pIBD[1] 
			J2 <- 0
		}
		#coe <- 0.5^(nmei-1); OMcoe <- 1. - coe
		if (all(pIBD<=1)) {
			#print(c(pIBD,J0,J1,J2))
			out <- array(c(
							J0*p2 + J1*p  + J2, J0*pq2 + J1*q         , J0*q2,
							J0*p2 + J1*p/2      , J0*pq2 + J1*(p+q)/2 + J2, J0*q2 + J1*q/2,
							J0*p2             , J0*pq2  + J1*p        , J0*q2  + J1*q + J2
					),			
					dim=c(3,3))
#			out <- array(c(
#							J0*p2*p2  + J1*p2  + J2, J0*p2*pq2  + J1*p*q         , J0*p2*q2,
#							J0*p2*pq2 + J1*p*q     , J0*pq2*pq2 + J1*(p2+q2) + J2, J0*q2*pq2 + J1*p*q,
#							J0*p2*q2               , J0*q2*pq2  + J1*p*q         , J0*q2*q2  + J1*q2 + J2
#					),			
#					dim=c(3,3))
#			out <- array(c(coe*p     + OMcoe*(p2), coe*q   + OMcoe*pq2,             OMcoe*(q2),
#							coe*0.5*p + OMcoe*(p2), coe*0.5 + OMcoe*pq2, coe*0.5*q + OMcoe*(q2),
#							OMcoe*(p2), coe*p   + OMcoe*pq2, coe*q     + OMcoe*(q2)),
#					dim=c(3,3))
#			out <- array(c(coe*p     + OMcoe*(p2), coe*q   + OMcoe*pq2,             OMcoe*(q2),
#							coe*0.5*p + OMcoe*(p2), coe*0.5 + OMcoe*pq2, coe*0.5*q + OMcoe*(q2),
#							OMcoe*(p2), coe*p   + OMcoe*pq2, coe*q     + OMcoe*(q2)),
#					dim=c(3,3))
#			out <- array(c(coe*p     + OMcoe*(p2+pq*F), coe*q   + OMcoe*pq2*(1-F),             OMcoe*(q2+pq*F),
#							coe*0.5*p + OMcoe*(p2+pq*F), coe*0.5 + OMcoe*pq2*(1-F), coe*0.5*q + OMcoe*(q2+pq*F),
#							OMcoe*(p2+pq*F), coe*p   + OMcoe*pq2*(1-F), coe*q     + OMcoe*(q2+pq*F)),
#					dim=c(3,3))
		} else {
			out <- array(c(1, 0, 0,
							0, 1, 0,
							0, 0, 1),
					dim=c(3,3))
		}
		return(out)
	}
	if (length(q)>1) {
		q <- array(q)
		out <- sapply(q,FUN=getTT,pIBD=.5^(nmeioses-1),simplify="array")
	} else {
		out <- getTT(q,pIBD=.5^(nmeioses-1))
	}
	return(out)
}
