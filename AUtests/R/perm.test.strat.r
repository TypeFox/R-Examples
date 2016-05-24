#' Stratified permutation testing
#'
#' Calculates permutation p-values for testing independence in 2x2 case-control tables, while adjusting for categorical covariates. 
#' Inputs are given as a vector of counts in each strata defined by the covariate(s).  
#' @param m0list Number of control subjects in each strata
#' @param m1list Number of case subjects in each strata
#' @param r0list Number of control subjects exposed in each strata
#' @param r1list Number of case subjects exposed in each strata
#' @return A permutation p-value, computed under the likelihood ratio test.
#' @examples
#' perm.test.strat(c(7000, 1000), c(11000, 1000), c(50, 30), c(70, 40))

perm.test.strat = function(m0list, m1list, r0list, r1list)
{
	q = length(m0list) # number of strata
		
	# For each strata:
	# - make data frame (dd) of all possible (permuted) outcomes
	# - calculate (log) probability of each outcome being observed
	# - calculate LR stat for each outcome
	dd.list = lapply(1:q, function(i)
	{
		dd = data.frame(r0x=0:(r0list[i]+r1list[i]))
		dd$r1x = r0list[i] + r1list[i] - dd$r0x

		dd$prob = dhyper(dd$r1x, r0list[i]+r1list[i], m0list[i]+m1list[i]-r0list[i]-r1list[i], m1list[i], log=TRUE)

		dd$p0x = with(dd, r0x/m0list[i])
		dd$p1x = with(dd, r1x/m1list[i])
		dd$pLx = with(dd, (r0x+r1x)/(m0list[i]+m1list[i]))

		dd$llik.nullx = with(dd, dbinom(r0x, m0list[i], pLx, log=T) + dbinom(r1x, m1list[i], pLx, log=T))
		dd$llik.altx = with(dd, dbinom(r0x, m0list[i], p0x, log=T) + dbinom(r1x, m1list[i], p1x, log=T))
		dd$llrx = with(dd, llik.altx - llik.nullx)
		
		dd
	})

	lengths = sapply(1:q, function(i) dim(dd.list[[i]])[1]) # number of outcomes in each strata

	# all combinations of outcomes across all strata (dim (l1*l2*...*lq) x q )
	dd2 = expand.grid(lapply(1:q, function(i) 1:lengths[i])) 
	# compute overall probabilities and LR test statistics
	dd2$prob = exp(rowSums( sapply(1:q,  function(i){dd.list[[i]]$prob[dd2[,i]]} ) ))
	dd2$llrx = rowSums( sapply(1:q,  function(i){dd.list[[i]]$llrx[dd2[,i]]} ) )

	# Compare LR stats to observed data
	matchrows = sapply(1:q, function(i) which(with(dd.list[[i]], (r0x==r0list[i]) & (r1x==r1list[i])))) # matches in dd.list (q matches)
	matchrow = which(apply(dd2[,1:q], 1, function(indexes){ all(indexes==matchrows) } )) # match in dd2 (1 match)

	# Observed LR stat
	llr = dd2[matchrow, "llrx"]
  p.obs = dd2[matchrow, "prob"] 
  dd2 = dd2[-matchrow,]
	
	# Permutation p-value: sum probabilities for all outcome combinations with LR stats >= observed LR stat
	dd.lrt = p.obs + sum(dd2[dd2$llrx >= llr,]$prob)
	c(lrt.p=dd.lrt)
}
