print.hogg.test <-
function(x, digits = max(5, .Options$digits - 2), ...) {
	cat("\nHogg\'s Adaptive Test\n")
	res<-c(x$statistic,x$p.value)
	names(res)<-c('Statistic','p.value')
	print(format(res,digits=digits),quote=FALSE)
	cat("\nScores Selected: ", x$scores, "\n")
}
