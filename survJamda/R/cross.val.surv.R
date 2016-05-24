cross.val.surv <-
function (x, y, censor, ngroup, iter, method, zscore,gn.nb,gn.nb.display,plot.roc)
{
   call <- match.call()
    x <- as.matrix(x)

    n <- length(y)

    if (ngroup < 2) {
        stop("ngroup should be greater than or equal to 2")
    }
    if (ngroup > n) {
        stop("ngroup should be less than or equal to the number of observations")
    }
    if (ngroup == n) {
        groups <- sample(1:n)
        leave.out <- 1
    }
    if (ngroup < n) 
	groups = groups.cv (n, ngroup,censor) 

     	all.group = NULL
   	pred.fit = vector (length = ngroup)
    	sign = NULL	
	all.fp = NULL
	all.tp = NULL
	options(warn=-1) 
   for (j in 1:ngroup) {
	if (zscore){
		x[-groups[[j]], ] = scale (t(scale(t(x[-groups[[j]], ]))))
		x[groups[[j]], ] = scale (t(scale(t(x[groups[[j]], ]))))
	}

	my.func <- featureselection
	p.list<- do.call(my.func, list(x[-groups[[j]],], y[-groups[[j]]],censor[-groups[[j]]], method, gn.nb))

		cox.coef = cal.cox.coef (x[-groups[[j]],], y[-groups[[j]]],censor[-groups[[j]]])

	if (gn.nb.display >= 1){
		if (is.matrix(x)){
			cat("Selected genes:\n")
			cat(sort(colnames(x)[p.list][1:gn.nb.display]), "\n")
		}
	}

	lp.train = cox.coef[p.list]%*%t(x[-groups[[j]],p.list])
	lp.train = as.vector(lp.train)

	if (is.vector(x[groups[[j]],p.list]) && length(x[groups[[j]],p.list]) == length(cox.coef[p.list]))
		m = x[groups[[j]],p.list]
	else
		m = t(x[groups[[j]],p.list])

	lp = cox.coef[p.list]%*%m
	lp = as.vector(lp)

	predict.time =  mean(y[groups[[j]]][censor[groups[[j]]]==1])

	roc.fit =survivalROC (Stime = y[groups[[j]]], status = censor[groups[[j]]], marker =lp, predict.time = predict.time, span = 0.25*NROW(x[groups[[j]],])^(-0.20))
	pred.fit[j] = roc.fit$AUC

	all.fp = rbind(all.fp,roc.fit$FP)
	all.tp = rbind(all.tp,roc.fit$TP)

	sgn = ifelse (lp < median(lp.train),0, 1)
	sgn = as.vector(sgn)
	sign = c(sign, sgn)
	all.group = c(all.group,groups[[j]])
    }

	if (plot.roc) lines(mean(as.data.frame(all.fp)),mean(as.data.frame(all.tp)), lty = 3)

	cox.hr = coxph(Surv(as.vector (y[all.group]),as.vector (censor[all.group]))~sign)
	cat ("Iteration", iter, "\t",sprintf("%.2f",mean(pred.fit)), "\t", sprintf("%.2f",summary (cox.hr)$coefficients[,"exp(coef)"]), "(", sprintf("%.2f",summary (cox.hr)$conf.int[,"lower .95"]), "-", sprintf("%.2f",summary (cox.hr)$conf.int[,"upper .95"]),")\t", summary (cox.hr)$coefficients[,"Pr(>|z|)"],"\n", sep = "")

	val = c(mean(pred.fit), summary (cox.hr)$coefficients[,"exp(coef)"])
    return(val)
}
