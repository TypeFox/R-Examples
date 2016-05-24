plotROC <-
function(test.ind,all.surv,all.censor, lp, file.name,col, normalization, time.dep)
{
	roc.fit = NULL
	surv = all.surv[test.ind]
	censor = all.censor[test.ind]

	file.name = detFileName(file.name)

	excl.ind = excl.samples(test.ind,surv,censor)
	
	if (is.null(excl.ind)){
		ind = test.ind
		lp.ind = 1:length(test.ind)
	}
	else{
		first.ind = test.ind[1]-1
		ind = setdiff(test.ind,first.ind+excl.ind)
		lp.ind = ind-first.ind
	}

	if (time.dep)
		plot.time.dep (all.surv[ind],all.censor[ind],lp[lp.ind],test.ind, file.name, col)
	else
		plot.roc.curves (all.surv[ind],all.censor[ind],lp[lp.ind],test.ind, file.name, col, normalization)
}

