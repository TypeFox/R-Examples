calPerformance.auc.plot <-
function(lst, train.ind, test.ind, file.name,col, method, normalization, time.dep)
{
	train = lst$mat[train.ind,]
	test = lst$mat[test.ind,]
	options(warn = -1)

        res = featureselection (train, lst$phyno$surv[train.ind], lst$phyno$censor[train.ind], method, 100)

	my.func <- featureselection
        list.p<- do.call(my.func, list(train, lst$phyno$surv[train.ind], lst$phyno$censor[train.ind], method, 100))

   	cox.coef = cal.cox.coef (train, lst$phyno$surv[train.ind], lst$phyno$censor[train.ind])

	lp.train = cox.coef[list.p]%*%t(train[,list.p])
	lp = cox.coef[list.p]%*%t(test[,list.p])

	plotROC(test.ind,lst$phyno$surv,lst$phyno$censor, lp, file.name,col, normalization, time.dep)

}

