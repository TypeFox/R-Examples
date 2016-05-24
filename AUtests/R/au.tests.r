# AU p-value
# * Input: single dataset (m0 = # control, m1 = # case, r0 = # control with allele, r1 = # case with allele)
# * Compute test statistic, T
# * Given n = m0+m1, p = (r0+r1)/(m0+m1), get all plausible datasets (m0, m1, r0x, r1x), assign P(r0x,r1x) = P(r0x)P(r1x) (null)
# * Compute test statistics Tx for all plausible datasets
# * Output: sum probabilities of datasets where abs(Tx) >= abs(T)

#' AU testing
#'
#' Calculates approximate unconditional p-values for testing independence in 2x2 case-control tables.
#' @param m0 Number of control subjects
#' @param m1 Number of case subjects
#' @param r0 Number of control subjects exposed
#' @param r1 Number of case subjects exposed
#' @param lowthresh A threshold for probabilities below to be considered as zero. Defaults to 1e-12.
#' @return A vector of AU p-values, computed under score, likelihood ratio, and Wald tests.
#' @examples
#' au.tests(15000, 5000, 30, 25)
#' au.tests(10000, 10000, 30, 25)

au.tests = function(m0, m1, r0, r1, lowthresh=1E-12) 
{	
 	if (r0 == 0 & r1 == 0)
	{
		return(c(score.p = 1, lr.p = 1, wald.p = 1, wald0.p = 1))
	}
	
 	if (r0 == m0 & r1 == m1)
	{
		return(c(score.p = 1, lr.p = 1, wald.p = 1, wald0.p = 1))
	}
		
	if (is.na(m0 + m1 + r0 + r1))
	{
		return(c(score.p = NA, lr.p = NA, wald.p = NA, wald0.p = NA))
	}

 	p = (r0+r1)/(m0+m1) # observed p
	
	# Score test
 	ybar = m1/(m0+m1)
 	t = r1*(1-ybar) - r0*ybar
 	sd.t = sqrt((1-ybar)^2*m1*p*(1-p) + ybar^2*m0*p*(1-p))

	# Likelihood ratio test	
 	p0 = r0/m0
 	p1 = r1/m1
 	pL = (r0+r1)/(m0+m1)

 	llik.null = dbinom(r0, m0, pL, log=T) + dbinom(r1, m1, pL, log=T)
 	llik.alt  = dbinom(r0, m0, p0, log=T) + dbinom(r1, m1, p1, log=T)
 	llr = llik.alt - llik.null

	# Wald test (with regularization)
 	reg = 0.5
 	betahat = log( (r1+reg)/(m1-r1+reg)/((r0+reg)/(m0-r0+reg)) )
 	sehat = sqrt(1/(r0+reg) + 1/(r1+reg) + 1/(m0-r0+reg) + 1/(m1-r1+reg))
 	waldT = betahat/sehat

	# Wald test (no regularization)
 	reg0 = 0
 	betahat0 = log( (r1+reg0)/(m1-r1+reg0)/((r0+reg0)/(m0-r0+reg0)) )
 	sehat0 = sqrt(1/(r0+reg0) + 1/(r1+reg0) + 1/(m0-r0+reg0) + 1/(m1-r1+reg0))
 	waldT0 = betahat0/sehat0
	
	# Approximate unconditional	p-value

 	hicount = qbinom(lowthresh, m0+m1, p, lower.tail=F)

 	dd = expand.grid(r0x=0:hicount, r1x=0:hicount)
	dd = dd[-1,]
 	dd$prob = dbinom(dd$r0x, m0, p)*dbinom(dd$r1x, m1, p)
		
 	dd$px = with(dd, (r0x+r1x+1)/(m0+m1+2) )
 	dd$tx = with(dd, r1x*(1-ybar) - r0x*ybar )
 	dd$sd.tx = with(dd, sqrt((1-ybar)^2*m1*px*(1-px) + ybar^2*m0*px*(1-px)))

 	dd$p0x = with(dd, r0x/m0)
 	dd$p1x = with(dd, r1x/m1)
 	dd$pLx = with(dd, (r0x+r1x)/(m0+m1))

 	dd$llik.nullx = with(dd, dbinom(r0x, m0, pLx, log=T) + dbinom(r1x, m1, pLx, log=T))
 	dd$llik.altx  = with(dd, dbinom(r0x, m0, p0x, log=T) + dbinom(r1x, m1, p1x, log=T))
 	dd$llrx = with(dd, llik.altx - llik.nullx)
	
	dd$betahatw = with(dd, log(r1x+ reg) - log(m1-r1x+reg) - log(r0x+reg) + log(m0-r0x+reg) )
 	dd$sehatw = with(dd, sqrt(1/(r0x+reg) + 1/(r1x+reg) + 1/(m0-r0x+reg) + 1/(m1-r1x+reg)) )
 	dd$waldTw = with(dd, betahatw/sehatw)
	
	dd$betahatw0 = with(dd, log(r1x+ reg0) - log(m1-r1x+reg0) - log(r0x+reg0) + log(m0-r0x+reg0) )
 	dd$sehatw0 = with(dd, sqrt(1/(r0x+reg0) + 1/(r1x+reg0) + 1/(m0-r0x+reg0) + 1/(m1-r1x+reg0)) )
 	dd$waldTw0 = with(dd, betahatw0/sehatw0)
	infrows = which( with(dd, abs(betahatw0) == Inf) )
  p.inf = sum(dd[infrows, "prob"])

	matchrow = which( with(dd, r0x==r0 & r1x==r1) )
  p.obs = dd[matchrow, "prob"]
  dd = dd[-matchrow,]
  
	
 	dd.score = p.obs + sum(dd[abs(dd$tx/dd$sd.tx) >= abs(t/sd.t),]$prob)
	dd.lrt   = p.obs + sum(dd[dd$llrx >= llr,]$prob)
 	dd.wald  = p.obs + sum(dd[abs(dd$waldTw) >= abs(waldT),]$prob)
  dd = dd[-which( with(dd, abs(betahatw0) == Inf) ),]
	dd.wald0 = p.obs + sum(dd[abs(dd$waldTw0) >= abs(waldT0),]$prob) + p.inf
	
  c(score.p = dd.score, lr.p = dd.lrt, wald.p=dd.wald, wald0.p = dd.wald0)
}

