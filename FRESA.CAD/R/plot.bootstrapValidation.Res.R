#' @method plot bootstrapValidation_Res
plot.bootstrapValidation_Res <-
function(x,xlab = "Years", ylab="Survival",...) 
{


	par(mfrow=c(1,1),pty='m')
	classlen=length(class(x$boot.model))
	
	cobj <- substr(class(x$boot.model)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			varsList <- as.list(attr(terms(x$boot.model),"variables"))
			form = paste(varsList[2]," ~ strata(x$boot.model$linear.predictors >0)")
			cat(form)
			plot(survival::survfit(formula(form),data = x$data),xlab=xlab,ylab=ylab,col=c("blue","red"),conf.int=TRUE,lty = 2:3,...)
			legend(0.1, 0.3, c("Low Risk", "High Risk"), col=c("blue","red"), lty = 2:3)
			title("Kaplan-Meier Curve");
		},
		{
#			plot(x$boot.model,...)
		}
	)
	pROC::roc(x$testOutcome,x$testPrediction,col="red",auc=TRUE,print.auc=TRUE,plot=TRUE,smooth=FALSE)
	par(new=TRUE)
	pROC::roc( x$outcome, x$boot.model$linear.predictors,plot=TRUE,ci=TRUE,auc=TRUE,of='se',specificities=c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),boot.n=200,smooth=FALSE);        
}
