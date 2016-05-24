


summary.ClaimsRules <-
function(object, ...)
{
	#assignInNamespace("summary.ClaimsRules", summary.ClaimsRules, ns = asNamespace("base"))
	x<-object
	cat("\n")
	cat("Claims of the Agents","\n")
	cat("\n")
	cl<-nrow(x)-1
	M1<-x[1:cl,1]
	print(M1)
	cat("\n")
	cat("Assignments according to the following rules","\n")
	cat("\n")
	M2<-x[1:cl,-1]
	print(M2)
	cat("\n")
	cat("Inequality Analysis among rules (Gini Index)","\n")
	cat("\n")
	M3<-x[cl+1,-1]
	print(M3)
}
