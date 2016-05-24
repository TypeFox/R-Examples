"sackin.test" <-
function(tree,model="yule",alternative="less",n.mc=500) {
	

	#catalan <- function(n) {
	# 	if (n<2) return(1)
	#	s<-seq(1,n,by=2)
	#	res<-prod(s)
	#	return(res)
	#}
	
	
	sackin3<-function(n,model,p=1/3) {
		if (n>2) {
			l <- switch( model,
			"biased" = sample(n-1,1,.5*dbinom(0:(n-2),size=(n-2),prob=p)+.5*dbinom(0:(n-2),size=(n-2),prob=(1-p))),
			"yule" = sample(n-1,1),
			"aldous" = sample(n-1,1, prob = (1/(1:(n-1))/((n-1):1)) ) )
			#"pda" = sample(n-1,1, prob=(sapply(1:(n-1), FUN=function(x) { 1/2*choose(n,x)*catalan(x)*catalan(n-x)/catalan(n) })))
			return( sackin3(l,model,p) + sackin3(n-l,model,p) + n ) 
		}
		else return(0)
	}
	
	
	
	
	
	#defender 
	if (class(tree)!='treeshape') {
		stop("invalid arguments")
	}
	
	if ((model != "yule") && (model!="pda")) {
		stop("Argument 'model' is incorrect. Should be a string 'pda' or 'yule'")
	}
	
	tip.number.mc<-nrow(tree$merge)+1
	if (model=="pda") {
		cat("Computing Monte Carlo estimates...")
		trees.mc<-rtreeshape(n=n.mc,tip.number=tip.number.mc,model=model)
		lind.mc<-sapply(trees.mc,FUN=sackin,norm=NULL)
	}
	else {
	
		lind.mc<-NULL
		for (i in 1:n.mc) {
			lind.mc<-c(lind.mc,sackin3(tip.number.mc,model))
		}
	}
	res<-ecdf(lind.mc) 
	cat("\n\n")

	#m<-mean(lind.mc)
	
	#resbis<-quantile(lind.mc,probs=seq(0,1,length=101))
	
	
	stat<-sackin(tree,norm=NULL)
	stat_norm<-sackin(tree,norm=model)
	cat("	Test of the ")
	cat(model)
	cat(" hypothesis using the Sackin index ","\n")
	cat("Statistic = ")
	cat(stat, "\n")
	cat("Standardized Statistic = ")
	cat(stat_norm,"\n")
	cat("p-value = ")
	if (alternative == "less") {
		p.value<-res(stat)
		cat(p.value,"\n")
		cat("alternative hypothesis: the tree is more balanced than predicted by the ")
		cat(model)
		cat(" model","\n") 
	}
	else if (alternative == "greater") {
		p.value<-1-res(stat)
		cat(p.value,"\n")	
		cat("alternative hypothesis: the tree is less balanced than predicted by the ")
		cat(model)
		cat(" model","\n") 
	}
	else {
		cat("\n")
		stop("Argument 'alternative' is incorrect. Should be a string 'less' or 'greater'")
	}
	cat("\n")
	cat("Note : the p-value was computed using a Monte-Carlo method")
	cat("\n")
	
	res<-list(model=model,statistic = stat,p.value = p.value, alternative = alternative)
	
}

