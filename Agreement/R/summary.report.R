summary.report <-
function(object,...){
	cat("Agreement Statistics\n");
	summary(object$Estimate,object$Data$dec,...);
	cat("\n");
	cat(object$Data$conf);cat("% Confidence Limits\n");
	summary(object$Conf_Limit,object$Data$dec,...);
	cat("\n");
	cat("Allowance\n");
	summary(object$Allowance,object$Data$dec,...);
}

