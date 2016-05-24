#' Firth AU testing
#'
#' Calculates approximate unconditional Firth test p-value for testing independence in 2x2 case-control tables.
#' The Firth test requires more computational resources than the tests computed in the au.tests function.
#' @param m0 Number of control subjects
#' @param m1 Number of case subjects
#' @param r0 Number of control subjects exposed
#' @param r1 Number of case subjects exposed
#' @param lowthresh A threshold for probabilities below to be considered as zero. Defaults to 1e-12.
#' @return A single AU p-value, computed under the Firth test.
#' @examples
#' au.firth(15000, 5000, 1, 0)

au.firth = function(m0, m1, r0, r1, lowthresh=1E-12) 
{	
 	if (r0 == 0 & r1 == 0)
	{
		return(c(au.firth.p = 1))
	}
	
 	if (r0 == m0 & r1 == m1)
	{
		return(c(au.firth.p = 1))
	}
		
	if (is.na(m0 + m1 + r0 + r1))
	{
		return(c(au.firth.p = NA))
	}

 	p = (r0+r1)/(m0+m1) # observed p
	
	y = c(1,1,0,0)
 	x = c(1,0,1,0)
 	weights = c(r1, m1-r1, r0, m0-r0)
	
 	firth.t = sum(c(-2,2)*logistf(y~x, weights=weights)$loglik)

	# Approximate unconditional p-value
 	hicount = qbinom(lowthresh, m0+m1, p, lower.tail = F)
	
	dd = expand.grid(r0x=0:hicount, r1x=0:hicount)
	dd = dd[-1,]
	dd = subset(dd, dd$r0x + dd$r1x <= hicount)	
 	dd$prob = dbinom(dd$r0x, m0, p)*dbinom(dd$r1x, m1, p)
	
	dd$firth.tx = sapply(1:nrow(dd), function(a)
	{sum(c(-2,2)*logistf(c(1,1,0,0) ~ c(1,0,1,0), weights = c(dd$r1x[a], m1-dd$r1x[a], dd$r0x[a], m0-dd$r0x[a]))$loglik)})
		
	matchrow = which( with(dd, r0x==r0 & r1x==r1) )
  p.obs = dd[matchrow, "prob"]
  dd = dd[-matchrow,]
 	dd.firth = p.obs + sum(dd[dd$firth.tx >= firth.t,]$prob)
	
  c(au.firth.p = dd.firth)
}
