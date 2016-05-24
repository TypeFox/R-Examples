eval.subset <-
function (x, y, censor,iter, method, gn.nb, train.nb)
{
   call <- match.call()
    x <- as.matrix(x)
    n <- length(y)
  
	gn.lst = NULL
	 options(warn=-1) 

      	lst.samples = shuffle.samples(n, censor, train.nb)	
	train.ind = lst.samples$train.ind
	test.ind = lst.samples$test.ind
	
	my.func <- featureselection
	p.list<- do.call(my.func, list(x[train.ind,], y[train.ind],censor[train.ind], method, gn.nb))
	
	cox.coef = cal.cox.coef(x[train.ind,], y[train.ind],censor[train.ind])
  	
	lp.train = cox.coef[p.list]%*%t(x[train.ind,p.list])
	lp.train = as.vector(lp.train)

	if (is.vector(x[test.ind,p.list]) && length(x[test.ind,p.list]) == length(cox.coef$coef[p.list]))
		m = x[test.ind,p.list]
	else
		m = t(x[test.ind,p.list])

	lp = cox.coef[p.list]%*%m
	lp = as.vector(lp)
    
	roc.fit =survivalROC (Stime = as.vector(y[test.ind]), status = as.vector(censor[test.ind]), marker=lp, predict.time = mean(y[test.ind]), span = 0.25*NROW(x[test.ind,])^(-0.20))

	sgn = ifelse (lp < median(lp.train),0, 1)
	sgn = as.vector(sgn)
	
	cox.hr = coxph(Surv(as.vector (y[test.ind]),as.vector (censor[test.ind]))~sgn)
        if(summary (cox.hr)[[6]][2] < 500){
		cat ("Iteration", iter, "\t",sprintf("%.2f",mean(roc.fit$AUC)), "\t", sprintf("%.2f",summary (cox.hr)$coefficients[,"exp(coef)"]), "(", sprintf("%.2f",summary (cox.hr)$conf.int[,"lower .95"]), "-", sprintf("%.2f",summary (cox.hr)$conf.int[,"upper .95"]),")\t", summary (cox.hr)$coefficients[,"Pr(>|z|)"],"\n", sep = "")

		val = c(mean(roc.fit$AUC), summary (cox.hr)$coefficients[,"exp(coef)"])
	}
	else
		 val = c(0,0)

    	return(val)
}

