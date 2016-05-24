
combiopt <- function(g){
	ind1 <- c(1, 1, 2, 3)
	ind2 <- c(2, 4, 3, 4)
	excess <- sum(g)/2 - g[ind1]- g[ind2]
	excessopt <- max(excess)
	ind1 <- ind1[excess==excessopt][1]
	ind2 <- ind2[excess==excessopt][1]
	return(list(ind1=ind1, ind2=ind2, excessopt=excessopt))
	}

ULg <- function(U, L){
	u1 <- U[1]
	u2 <- U[2]
	u3 <- U[3]
	l1 <- L[1]
	l2 <- L[2]
	#the excess
	s <- l1 + l2 - u1 - u2 - u3 
	g1 <- l1 + l2 - u2 - u3
	g2 <- l1 + l2 - u1 - u3
	g3 <- l1 + l2 - u1 - u2
	g <- c(g1, g2, g3, l1, l2)
	return(list(g=g, excess=s))
}


Thomae <- function(U, L, lB, tol, maxiter, debug){
	V <- ULg(U,L)
	g <- V$g                             # the permuting variables in Thomae's theorem
	s <- V$excess                        # the excess corresponding to the input permutation
	if (s <= 0) return(list(G1=NA))      # negative excess
	M <- combiopt(g)                     # optimal combination of Thomae's arguments
	excessopt <- M$excessopt
	Lopt <- g[c(M$ind1,M$ind2)]	         # lower parameters corresponding to the maximum excess
	Uopt <- g[-c(M$ind1,M$ind2)] - excessopt
	Gg <- genhypergeo_series(Uopt,Lopt,1, tol=tol, maxiter=maxiter, debug=debug)
	F32 <- Gg[[1]]
	out <- NULL
	if (debug) out <- Gg[[2]]
# Equivalence factors are shape1 product - ratio of gammas.
# First compute at the logarithmic scale to avoid numerical problems with large values
# The log hypergeometric 3F2(U,L;1) is then given by:
	logG1 <- sum(lgamma(L)-lgamma(Lopt))+ lgamma(s) -lgamma(excessopt) + log(F32)
# The term in the Gini expression is thus
	logG1 <- logG1 + lB 
	G1 <- exp(logG1)
	return(list(G1=G1, F32=F32, Uopt=Uopt, Lopt=Lopt, out=out))
}


gb2.gini <- function(shape1, shape2, shape3, tol=1e-08, maxiter=10000, debug=FALSE){
  if (shape1 < 0 | shape2 < 0 | shape3 < 0) {print("Warning: negative parameter", quote=FALSE); return(list(G1=NA, G2=NA, B=NA, Gini=NA))}
	excess <- shape3-1/shape1
	if (excess <= 0) {print("Warning: non-positive excess; expectation does not exist", quote=FALSE); return(list(G1=NA, G2=NA, B=NA, Gini=NA))}
	if (excess < 1e-10) {print("Warning: excess less than 1e-10; limiting value of 1 forced for Gini", quote=FALSE)
			     return(list(G1=NA, G2=NA, B=NA, Gini=1))}
	U <- c(1, shape2 + shape3, 2*shape2 + 1/shape1)
	L1 <- c(shape2 + 1, 2*(shape2 + shape3))
	L2 <- c(shape2 + 1 + 1/shape1, 2*(shape2 + shape3))
	lB <- lbeta(2*shape3-1/shape1, 2*shape2+1/shape1) - lbeta(shape2, shape3)-lbeta(shape2 + 1/shape1, shape3 - 1/shape1)
	T1 <- Thomae(U, L1, lB, tol, maxiter, debug)
	T2 <- Thomae(U, L2, lB, tol, maxiter, debug)
	G1 <- T1$G1
	G2 <- T2$G1
	Gini <- G1/shape2 - G2/(shape2+1/shape1)
	if(is.infinite(G1)) {print("Warning: overflow occured; limiting value of 1 forced for Gini", quote=FALSE)
                             return(list(G1=NA, G2=NA, B=NA, Gini=1))}
	if(!is.na(Gini) & Gini > 1+1e-10)
	{print("Warning! Gini estimate > 1:", quote=FALSE); print(Gini); print("Gini forced to 1", quote=FALSE); Gini <-1}
	return(list(G1=G1, G2=G2, F321 = T1$F32, F322= T2$F32, lB=lB, Gini=Gini, U=U, L1=L1, L2=L2, Uopt1=T1$Uopt, Lopt1=T1$Lopt, T1out=T1$out, T2out=T2$out))
}
