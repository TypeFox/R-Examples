print.rank.test<-function(x,...) {
	with(x,cat("statistic = ",statistic,", p-value = ", p.value, "\n"))
	if(with(x,exists('conf.int'))) {
		cat(100*x$conf.level," percent confidence interval:\n")
		cat(x$conf.int, "\n")
	}
	if(with(x,exists('estimate'))) {
		cat("Estimate:", x$estimate,"\n")	
	}
}

print.fkk.test<-function(x,...) {
	coef<-cbind(x$estimate,x$conf.int)
	colnames(coef)<-c('estimate','ci.lower','ci.upper')
	cat('Table of estimates and ', 100*x$conf.level," percent confidence intervals:\n")
	print(coef)
	with(x,cat("\n","Test statistic = ",statistic," p-value = ", p.value, "\n"))
}
