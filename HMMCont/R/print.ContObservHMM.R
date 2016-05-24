print.ContObservHMM <-
function(x, ...)
{
	cat("\nThe number of Baum-Welch iterations: ")
	cat((nrow(x$Parameters))-1)
	cat("\nThe parameters accumulated so far: \n")
	par<-round(x$Parameters, 2)
	print(par)
	cat("\nThe results accumulated so far: \n")
	res1<-(x$Results[-(nrow(x$Parameters)),3])
	res2<-(x$Results[-(nrow(x$Parameters)),4])
	res3<-(x$Results[-(nrow(x$Parameters)),5])
	res<-cbind(res1, res2, res3)
	colnames(res)<-c("P", "AIC", "SBIC")
	print(res)
	yesno<-ifelse(((x$Viterbi[1,1]) == 0), "not yet ", "already ") 
	cat("\nThe Viterbi algorithm was ")
	cat(yesno)
	cat("executed\n\n")
	
}
