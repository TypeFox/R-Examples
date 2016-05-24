"likelihood.test" <-
function(tree,model="yule",alternative="two.sided") {
	
	if (class(tree)[[1]]!='treeshape') {
		stop("invalid arguments")
	}

	#number of internal nodes
	n<-nrow(tree$merge)
	if (n<4) {
		stop("This test cannot be computed for trees with less than 4 leaves (negative variance)")
	}
	if (model=="yule") {
		stat<-(shape.statistic(tree, norm="yule"))
		cat("Test of the Yule hypothesis: \n")
		cat("statistic = ")
		cat(stat,"\n")
		if (alternative=="two.sided") {
			cat("p.value = ")
			p.value<-2*(1-pnorm(abs(stat)))
			cat(p.value,"\n") 
			cat("alternative hypothesis: the tree does not fit the Yule model")
			cat("\n")
		}
		else if (alternative=="less") {
			cat("p.value = ")
			p.value<-pnorm(stat)
			cat(p.value,"\n")
			cat("alternative hypothesis: the tree is more balanced than predicted by the Yule model")
			cat("\n")
		}
		else if (alternative=="greater") {
			cat("p.value = ")
			p.value<-1-pnorm(stat)
			cat(p.value,"\n")
			cat("alternative hypothesis: the tree is less balanced than predicted by the Yule model")
			cat("\n")
		}
		else {
			stop("alternative hypothesis invalid")
		}
	}
	else if (model=="pda") {
		stat<-(shape.statistic(tree, norm="pda"))
		cat("Test of the PDA hypothesis: \n")
		cat("statistic = ")
		cat(stat,"\n")
		if (alternative=="two.sided") {
			cat("p.value = ")
			p.value<-2*(1-pnorm(abs(stat)))
			cat(p.value,"\n") 
			cat("alternative hypothesis: the tree does not fit the PDA model")
			cat("\n")
		}
		else if (alternative=="less") {
			cat("p.value = ")
			p.value<-pnorm(stat)
			cat(p.value,"\n")
			cat("alternative hypothesis: the tree is more balanced than predicted by the PDA model")
			cat("\n")
		}
		else if (alternative=="greater") {
			cat("p.value = ")
			p.value<-1-pnorm(stat)
			cat(p.value,"\n")
			cat("alternative hypothesis: the tree is less balanced than predicted by the PDA model")
			cat("\n")
		}
		else {
			stop("alternative hypothesis invalid")
		}
	}
	else {
		stop("model invalid. Should be 'yule' or 'pda'")
	}
	cat("\n")
	cat("Note: the p.value was computed according to a normal approximation\n")
	res<-list(model=model,statistic=stat,p.value=p.value,alternative=alternative)
	
}

