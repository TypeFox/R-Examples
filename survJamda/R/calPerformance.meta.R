calPerformance.meta <-
function(common.gene,zstat, i, j, geno.files, surv.data, method)
{
	comb = inv.normal (i, zstat)

	total.p <- 2*pnorm(abs(comb),lower.tail=FALSE)
	p.list <- p.adjust(total.p,method=method)
	if (method == "none")
		p.list = order(p.list)[1:100]
	else
		p.list = p.list<= .05	

	names(comb) = common.gene
	cm.gn = common.gene[p.list]
	cm.gn = intersect(cm.gn, colnames(get(geno.files[j])))

	if (sum(p.list) > 0){
                lst.pheno = comb.surv.censor(geno.files,c(j), surv.data)
                surv = lst.pheno$surv
		censor = lst.pheno$censor

		test= get(geno.files[j])[!is.na(surv),]
		censor= censor[!is.na(surv)]
		surv= surv[!is.na(surv)]

		lp = comb[cm.gn]%*%t(test[,cm.gn])
		lp = as.vector(lp)

		sgn = ifelse(lp < median(lp), 0,1)
		sgn = as.vector(sgn)

		cox.hr = coxph(Surv(surv,censor)~sgn)

		roc.fit =survivalROC (Stime = surv, status = censor, marker=lp, predict.time = mean(surv), span = 0.25*NROW(test)^(-0.20))
		cat ("AUC\tHR(CI)\t\t\tP-val\n")
		cat (sprintf("%.2f",roc.fit$AUC), "\t", sprintf("%.2f",summary (cox.hr)$coefficients[,"exp(coef)"]), "(", sprintf("%.2f",summary (cox.hr)$conf.int[,"lower .95"]), "-", sprintf("%.2f",summary (cox.hr)$conf.int[,"upper .95"]),")\t\tp=", summary (cox.hr)$coefficients[,"Pr(>|z|)"],"\n",sep = "")
	}
}

