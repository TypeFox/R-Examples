calPerformance.single.indep <-
function(lst1, lst2, method,gn.nb,perf.eval)
{
	train.vec = lst1$phyno
	test.vec = lst2$phyno
	options(warn = -1)

	#list.p = fs (lst1$mat, train.vec$surv, train.vec$censor, method, gn.nb)
	
	my.func <- featureselection
	list.p<- do.call(my.func, list(lst1$mat, train.vec$surv, train.vec$censor, method, gn.nb))

	cox.coef = cal.cox.coef (lst1$mat, train.vec$surv, train.vec$censor)

	lp.train = cox.coef[list.p]%*%t(lst1$mat[,list.p])
	lp = cox.coef[list.p]%*%t(lst2$mat[,list.p])
	
	lp.train = as.vector(lp.train)
	lp = as.vector(lp)

	if(perf.eval == "auc"){
		roc.fit =survivalROC (Stime = test.vec$surv, status = test.vec$censor, marker=lp, predict.time = mean(test.vec$surv), span = 0.25*NROW(lst2$mat)^(-0.20))

		sgn = ifelse (lp < median(lp.train),0, 1)
		sgn = as.vector(sgn)

		cox.hr = coxph(Surv(test.vec$surv, test.vec$censor)~sgn)

		cat ("AUC\tHR(CI)\t\t\tP-val\n")
    
		cat (sprintf("%.2f",roc.fit$AUC), "\t", sprintf("%.2f",summary (cox.hr)$coefficients[,"exp(coef)"]), "(", sprintf("%.2f",summary (cox.hr)$conf.int[,"lower .95"]), "-", sprintf("%.2f",summary (cox.hr)$conf.int[,"upper .95"]),")\t\tp=", summary (cox.hr)$coefficients[,"Pr(>|z|)"],"\n",sep = "")
	}
	else
		if (perf.eval == "cindex"){
			res = concordance.index(x=lp,surv.time = test.vec$surv, surv.event =test.vec$censor)
			print(res)
		}
		else{
			dt.tr = data.frame ("time"=train.vec$surv, "event"=train.vec$censor, "score" = lp.train)
			dt.ts = data.frame ("time"=test.vec$surv, "event"=test.vec$censor, "score" = lp)
			res = sbrier.score2proba(dt.tr,dt.ts)
			print(res)
		}
    
}

