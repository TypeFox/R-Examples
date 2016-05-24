#' Basic testing
#'
#' Calculates standard p-values for testing independence in 2x2 case-control tables.
#' @param m0 Number of control subjects
#' @param m1 Number of case subjects
#' @param r0 Number of control subjects exposed
#' @param r1 Number of case subjects exposed
#' @return A vector of p-values, computed under score, likelihood ratio, Wald, Firth, and Fisher's exact tests.
#' @examples
#' basic.tests(15000, 5000, 30, 25)

basic.tests = function(m0, m1, r0, r1) 
{	
 	if (r0 == 0 & r1 == 0)
	{
		return(c(score.p = 1, lr.p = 1, wald.p = 1, wald0.p = 1, firth.p = 1, fisher.p = 1))
	}

 	if (r0 == m0 & r1 == m1)
	{
		return(c(score.p = 1, lr.p = 1, wald.p = 1, wald0.p = 1, firth.p = 1, fisher.p = 1))
	}

	if (is.na(m0 + m1 + r0 + r1))
	{
		return(c(score.p = NA, lr.p = NA, wald.p = NA, wald0.p = NA, firth.p = NA, fisher.p = NA))
	}
		
	p = (r0+r1)/(m0+m1) # observed p
	
	# Score test
 	ybar = m1/(m0+m1)
 	t = r1*(1-ybar) - r0*ybar
 	sd.t = sqrt((1-ybar)^2*m1*p*(1-p) + ybar^2*m0*p*(1-p))
	score.t = abs(t/sd.t)
	score.pv = 2*(1-pnorm(score.t))

	# Likelihood ratio test	
 	p0 = r0/m0
 	p1 = r1/m1
 	pL = (r0+r1)/(m0+m1)
 	llik.null = dbinom(r0, m0, pL, log=T) + dbinom(r1, m1, pL, log=T)
 	llik.alt  = dbinom(r0, m0, p0, log=T) + dbinom(r1, m1, p1, log=T)
 	llr = llik.alt - llik.null
	lr.pv = 1-pchisq(2*llr, df = 1)

	# Wald test (with regularization)
 	reg = 0.5
 	betahat = log( (r1+reg)/(m1-r1+reg)/((r0+reg)/(m0-r0+reg)) )
 	sehat = sqrt(1/(r0+reg) + 1/(r1+reg) + 1/(m0-r0+reg) + 1/(m1-r1+reg))
 	waldT = betahat/sehat
	wald.pv = 2*(1-pnorm(abs(waldT)))

	# Wald test (no regularization)
	if (r0 == 0 | r1 == 0) wald0.pv = 1
	else
	{
		reg0 = 0
		betahat0 = log( (r1+reg0)/(m1-r1+reg0)/((r0+reg0)/(m0-r0+reg0)) )
		sehat0 = sqrt(1/(r0+reg0) + 1/(r1+reg0) + 1/(m0-r0+reg0) + 1/(m1-r1+reg0))
		waldT0 = betahat0/sehat0
		wald0.pv = 2*(1-pnorm(abs(waldT0)))
	}

	# Firth test
	m.firth = logistf(c(1,1,0,0) ~ c(1,0,1,0), weights=c(r1, m1-r1, r0, m0-r0))
	firth.t = sum(c(-2, 2)*m.firth$loglik)
	firth.pv = m.firth$prob[2]
	names(firth.pv) = NULL
	
	# Fisher's exact test
	table.obs = matrix(c(m0-r0, m1-r1, r0, r1), nrow=2, 
	dimnames=list(Disease = c("Control", "Case"), Gene = c(0, 1)))
	fisher.pv = fisher.test(table.obs)$p.value

  	c(score.p = score.pv, lr.p = lr.pv, wald.p = wald.pv, wald0.p = wald0.pv, firth.p = firth.pv, fisher.p = fisher.pv)
}
