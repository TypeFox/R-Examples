# performs overall proportion test; if it is significant at alpha=0.05, performs pairwise tests to see where there are significant differences between proportions
# n should not be very different between the different groups under test; this is intended for simple analysis of split tests with multiple outcomes
# result rows are in the same order as the input x and n
significance_analysis <-
function(x, n) {
	num_groups = length(x)

	p = x/n
	my_rank = rank(-p)

	pt = prop.test(x=x, n=n)

	lower = rep(NA, num_groups)
	upper = rep(NA, num_groups)
	significance = rep(NA, num_groups)
	best = rep(0, num_groups)

	b = best_binomial_bandit(x, n)

	if (pt$p.value < 0.05) {
		is_best = 1
		for (cur_rank in sort(unique(my_rank))) {
			cur_indices = which(my_rank==cur_rank)
			if (length(cur_indices) >= 1) {
				cur_index = max(cur_indices)
				ranks_above = my_rank[my_rank > cur_rank]
				if (length(ranks_above) > 0) {
					comparison_index = min(which(my_rank==min(ranks_above)))
					# compare to the next lower proportion
					pt = prop.test(x=x[c(cur_index, comparison_index)], n=n[c(cur_index, comparison_index)], conf.level = (1 - 0.05))
					significance[cur_indices] = pt$p.value
					lower[cur_indices] = pt$conf.int[1]
					upper[cur_indices] = pt$conf.int[2]
					best[cur_indices] = is_best
					if (pt$p.value < 0.05) {
						is_best = 0
					}
				}
			}
		}
	}

	return(data.frame(successes=x, totals=n, estimated_proportion=p, lower=lower, upper=upper, significance=significance, rank=my_rank, best=best, p_best=b))
}
