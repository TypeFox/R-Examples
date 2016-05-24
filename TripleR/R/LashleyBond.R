# Compute LashleyBond standard errors for a single group
# @param saa, sbb, scc, sab, sccs = actor, partner, rel.ship-variance estimate, a-p-cov, rel-cov
# @param n = group size
compute_univariate_LB_SE2 <- function(saa, sbb, scc, sab, sccs, n) {
	
	# w is a general weighting factor
	w <- (n^2 - 3*n + 6) * (n^2 - 3*n + 4)

	sesaa2 <- (2*saa^2) / (n+1) + (2*(n^6 - 7*n^5 + 28*n^4 - 66*n^3 + 102*n^2 - 84*n + 32)* scc^2)/ (w*(n+1)*n^2*(n-2)^2) + (2*(n^3-n^2-2*n+16)*(n^2-2*n+2)*sccs^2) / (w*(n+1)*n^2*(n-2)^2) + (4*saa*((n-1) * scc + sccs)) / ((n+1)*n*(n-2)) + (4*(n^5-5*n^4+20*n^3-42*n^2+60*n-32)*scc*sccs) / (w*(n+1)*n^2*(n-2)^2)

	sesbb2 <- (2*sbb^2) / (n+1) + (2*(n^6 - 7*n^5 + 28*n^4 - 66*n^3 + 102*n^2 - 84*n + 32)* scc^2)/ (w*(n+1)*n^2*(n-2)^2) + (2*(n^3-n^2-2*n+16)*(n^2-2*n+2)*sccs^2) / (w*(n+1)*n^2*(n-2)^2) + (4*sbb*((n-1) * scc + sccs)) / ((n+1)*n*(n-2)) + (4*(n^5-5*n^4+20*n^3-42*n^2+60*n-32)*scc*sccs) / (w*(n+1)*n^2*(n-2)^2)

	sescc2 <- (2*(n^2-3*n+5)* ((scc^2 + sccs^2))) / w + (4*scc*sccs)/ w

	sesab2 <- ((n-3)*sab^2) / ((n+1)*(n-2)) + ((n^6 - 5*n^5 + 19*n^4 - 45*n^3 + 90*n^2 - 96*n + 64)* scc^2)/ (w*(n+1)*n^2*(n-2)^2) + ((n^6 - 7*n^5 + 31*n^4 - 83*n^3 + 150*n^2 - 144*n + 64)* sccs^2)/ (w*(n+1)*n^2*(n-2)^2) + ((n-1)*saa*sbb)/((n+1)*(n-2)) + ((n-1)*(saa+sbb)*((n-1)*scc+sccs))/((n+1)*n*(n-2)^2) + (2*(n-3)*sab*(scc+(n-1)*sccs))/((n+1)*n*(n-2)^2) + (4*(n^5-5*n^4+20*n^3-42*n^2+60*n-32)*scc*sccs) / (w*(n+1)*n^2*(n-2)^2)

	sesccs2 <- (2*(n^2-3*n+5)*((scc^2+sccs^2)))/w + (4*scc*sccs)/w
	
	return(c(sesaa2=sesaa2, sesbb2=sesbb2, sescc2=sescc2, NA, sesab2=sesab2, sesccs2=sesccs2))
}



# Compute LashleyBond standard errors for a single group
# @param saa, sbb, scc, sab, sccs = actor, partner, rel.ship-variance estimate, a-p-cov, rel-cov
# @param n = group size
compute_bivariate_LB_SE2 <- function(varComp.1, varComp.2, saf, sag, sbg, sbf, sch, schs, n) {

	coef1 <- n^11-14*n^10+89*n^9-342*n^8+872*n^7-1505*n^6+1698*n^5-1063*n^4+116*n^3+292*n^2-224*n+64
	coef2 <- n^10-12*n^9+66*n^8-227*n^7+534*n^6-857*n^5+883*n^4-416*n^3-148*n^2+224*n-64
	coef3 <- n^10-11*n^9+48*n^8-93*n^7-2*n^6+388*n^5-763*n^4+572*n^3+4*n^2-224*n+64
	coef4 <- n^11-16*n^10+117*n^9-520*n^8+1540*n^7-3083*n^6+3970*n^5-2689*n^4-4*n^3+1212*n^2-544*n-64
	coef5 <- n^10-11*n^9+42*n^8-33*n^7-258*n^6+976*n^5-1453*n^4+788*n^3+348*n^2-544*n-64
	coef6 <- 2*(n^10-14*n^9+92*n^8-383*n^7+1074*n^6-1963*n^5+2101*n^4-752*n^3-780*n^2+544*n+64)
	coef7 <- (n-3)*n*(n^6-9*n^5+35*n^4-75*n^3+76*n^2-12*n-48)

	suppressWarnings({
		sesaf2 <- ((n-1)*varComp.1[1]*varComp.2[1]+(n-3)*saf^2)/((n-2)*(n+1)) + ((n-1)^2*(varComp.1[1]*varComp.2[3]+varComp.1[3]*varComp.2[1])+(n-1)*(varComp.1[1]*varComp.2[6]+varComp.1[6]*varComp.2[1])+2*(n-3)*saf*((n-1)*sch+schs))/((n-2)^2*n*(n+1)) + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[6]+varComp.1[6]*varComp.2[3])+coef3*varComp.1[6]*varComp.2[6]+coef4*sch^2+coef5*schs^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7)

	
		sesbg2 <- ((n-1)*varComp.1[2]*varComp.2[2]+(n-3)*sbg^2)/((n-2)*(n+1)) + ((n-1)^2*(varComp.1[2]*varComp.2[3]+varComp.1[3]*varComp.2[2])+(n-1)*(varComp.1[2]*varComp.2[6]+varComp.1[6]*varComp.2[2])+2*(n-3)*sbg*((n-1)*sch+schs))/((n-2)^2*n*(n+1)) + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[6]+varComp.1[6]*varComp.2[3])+coef3*varComp.1[6]*varComp.2[6]+coef4*sch^2+coef5*schs^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7)
	
	sesch2 <- ((n^6-9*n^5+32*n^4-57*n^3+43*n^2+6*n-8)*(varComp.1[3]*varComp.2[3]+varComp.1[6]*varComp.2[6])+2*(n^4-6*n^3+3*n^2+18*n-8)*sch*schs)/coef7 + ((n^4-6*n^3+11*n^2-6*n+8)*(varComp.1[3]*varComp.2[6]+varComp.1[6]*varComp.2[3])+(n^6-9*n^5+28*n^4-33*n^3-9*n^2+54*n+8)*(sch^2 + schs^2))/coef7
		
	
	sesag2 <- ((n-1)*varComp.1[1]*varComp.2[2]+(n-3)*sag^2)/((n-2)*(n+1)) + ((n-1)^2*(varComp.1[1]*varComp.2[3]+varComp.1[3]*varComp.2[2])+(n-1)*(varComp.1[1]*varComp.2[6]+varComp.1[6]*varComp.2[2])+2*(n-3)*sag*((n-1)*schs+sch))/((n-2)^2*n*(n+1)) + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[6]+varComp.1[6]*varComp.2[3])+coef3*varComp.1[6]*varComp.2[6]+coef4*schs^2+coef5*sch^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7)

	
	sesbf2 <- ((n-1)*varComp.1[2]*varComp.2[1]+(n-3)*sbf^2)/((n-2)*(n+1)) + ((n-1)^2*(varComp.1[2]*varComp.2[3]+varComp.1[3]*varComp.2[1])+(n-1)*(varComp.1[2]*varComp.2[6]+varComp.1[6]*varComp.2[1])+2*(n-3)*sbf*((n-1)*schs+sch))/((n-2)^2*n*(n+1)) + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[6]+varComp.1[6]*varComp.2[3])+coef3*varComp.1[6]*varComp.2[6]+coef4*schs^2+coef5*sch^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7)


	seschs2 <- ((n^6-9*n^5+32*n^4-57*n^3+43*n^2+6*n-8)*(varComp.1[3]*varComp.2[3]+varComp.1[6]*varComp.2[6])+2*(n^4-6*n^3+3*n^2+18*n-8)*sch*schs)/coef7 + ((n^4-6*n^3+11*n^2-6*n+8)*(varComp.1[3]*varComp.2[6]+varComp.1[6]*varComp.2[3])+(n^6-9*n^5+28*n^4-33*n^3-9*n^2+54*n+8)*(sch^2 + schs^2))/coef7
	
	})
	
	return(c(sesaf2=sesaf2, sesbg2=sesbg2, sesag2=sesag2, sesbf2=sesbf2, sesch2=sesch2, seschs2=seschs2))
}