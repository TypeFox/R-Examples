calPerformance.merge.indep <-
function(lst, train.ind, test.ind, method,gn.nb,perf.eval, normalization)
{
	train = lst$mat[train.ind,]
	test = lst$mat[test.ind,]
	options(warn = -1)
	
	my.func <- featureselection
	list.p<- do.call(my.func, list(train, lst$phyno$surv[train.ind], lst$phyno$censor[train.ind], method, gn.nb))
       	cox.coef = cal.cox.coef (train, lst$phyno$surv[train.ind], lst$phyno$censor[train.ind])

	lp.train = cox.coef[list.p]%*%t(train[,list.p])
	lp = cox.coef[list.p]%*%t(test[,list.p])

	lp.train = as.vector(lp.train)
	lp = as.vector(lp)

	if(perf.eval == "auc"){
		roc.fit =survivalROC (Stime = lst$phyno$surv[test.ind], status =lst$phyno$censor[test.ind], marker=lp, predict.time = mean(lst$phyno$surv[test.ind]), span = 0.25*NROW(test)^(-0.20))

		sgn = ifelse (lp < median(lp.train),0, 1)
		sgn = as.vector(sgn)

		cox.hr = coxph(Surv(lst$phyno$surv[test.ind],lst$phyno$censor[test.ind])~sgn)

		cat ("AUC\tHR(CI)\t\tP-val\n")
		cat (sprintf("%.2f",roc.fit$AUC), "\t", sprintf("%.2f",summary (cox.hr)$coefficients[,"exp(coef)"]), "(", sprintf("%.2f",summary (cox.hr)$conf.int[,"lower .95"]), "-", sprintf("%.2f",summary (cox.hr)$conf.int[,"upper .95"]),")\t\tp=", summary (cox.hr)$coefficients[,"Pr(>|z|)"],"\n",sep = "")
	}
	else
		if (perf.eval == "cindex"){
			res = concordance.index(x=lp,surv.time = lst$phyno$surv[test.ind], surv.event =lst$phyno$censor[test.ind])
			print(res)
		}
		else{
			dt.tr = data.frame ("time"=lst$phyno$surv[train.ind], "event"=lst$phyno$censor[train.ind], "score" = lp.train)
			dt.ts = data.frame ("time"=lst$phyno$surv[test.ind], "event"=lst$phyno$censor[test.ind], "score" = lp)
			res = sbrier.score2proba(dt.tr,dt.ts)
			print(res)
		}
}

