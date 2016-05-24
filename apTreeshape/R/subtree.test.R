"subtree.test" <-
function(tree,size,alternative="two.sided") {
	#DEFENDER
	if (class(tree)!='treeshape') {
		stop("invalid arguments")
	}
	
	if (missing(size)) {
		size=2:(nrow(tree$merge))
	}
	#Mean
	m<-2/(size*(size+1))
	#Variance
	var<-2*(size-1)*(4*size^2-3*size-4)/(size*(2*size+1)*(2*size-1)*(size+1)^2)
	
	#frequency of substree of given size
	t<-spectrum.treeshape(tree)[nrow(tree$merge)-size+2]/(nrow(tree$merge)+1)
	
	#the following variable is asymptotically Gaussian for large tree sizes (mean=0 var=1)
	stat<-sqrt(nrow(tree$merge)+1)*(t-m)/sqrt(var)
	cat("     Test of the Yule hypothesis based on the number of subtrees of size ")
	cat(size,"\n")
	cat("Statistic = ")
	cat(t, "\n")
	cat("Standardized Statistic = ")
	cat(stat,"\n")
	cat("p-value = ")
	if (alternative == "two.sided") {
		p.value<-2*(1-pnorm(abs(stat)))
		cat(p.value)
		cat("\n")
		cat("alternative hypothesis: ")
		cat("The tree does not fit the Yule model","\n")
	}
	else if (alternative == "less") {
		p.value<-pnorm(stat)
		cat(p.value)
		cat("\n")
		cat("alternative hypothesis: ")
		cat("The tree is less balanced than predicted by the Yule model","\n")
	}
	else if (alternative == "greater") {
		p.value<-1-pnorm(stat)
		cat(p.value)
		cat("\n")
		cat("alternative hypothesis: ")
		cat("The tree is more balanced than predicted by the Yule model","\n")
	}
	else {
		cat("\n")
		stop("alternative argument invalid")
	}
	cat("\n")
	cat("Note : the p-value was computed according to a normal approximation.")
	cat("\n")
	res<-list(statistic = stat,p.value = p.value, alternative = alternative)
	
}

