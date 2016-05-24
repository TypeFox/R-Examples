getSummary.StratSel <-
function(obj,alpha=0.05, ...){
			if(missing(alpha)) alpha <- 0.05
			s <- summary.StratSel(obj)
			est <- s$coef[,1]
			se <- s$coef[,2]	
			stat <- s$coef[,3]
			p <- s$coef[,4]
			lwr <- est - qnorm((1-alpha/2))*se
			upr <- est + qnorm((1-alpha/2))*se
			N <- s$df
			coef <- cbind(est, se, stat, p, lwr, upr)
    		colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
		return(list(coef=coef, sumstat=c(N=N,AIC =s$AIC, AIC.c=s$AIC.c), call=obj$call))
}
