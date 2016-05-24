summary.StratSel <-
function(object, ...){
	se <- sqrt(diag(object$vcov))
	zval <- coef(object)/se
	
	TAB <- cbind(Estimate = coef(object),
                 StdErr = se,
                 z.value = zval,
                 p.value = 2*(1-pnorm(abs(zval))))
    colnames(TAB) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
   
    aic <- -round(2*object$logLik+2*length(object$coefficients),3)             
    aic.c <- -round(2*object$logLik+2*length(object$coefficients)+(2*length(object$coefficients)*(length(object$coefficients)+1))/(object$df-length(object$coefficients)-1),3)  
   		# aic.c based on Hurvich and Tsai 1989
   		# should be used when n/p<40
    res <- list(df=object$df,coefficients=TAB,ll=object$logLik,AIC=aic,AIC.c=aic.c, call=object$call, nits=object$nits, conv=object$conv)
    
    class(res) <- "summary.StratSel"
    res
}
