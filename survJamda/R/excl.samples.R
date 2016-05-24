excl.samples <-
function(test.ind,surv,censor)
{
	excl.ind = NULL
	for (w in (1:length(test.ind)))
		if (censor[which.min(surv)] == 0){
			excl.ind = c(excl.ind, which.min(surv))
			surv[which.min(surv)] = 1000
		}
		else
			break
	return (excl.ind)

}

