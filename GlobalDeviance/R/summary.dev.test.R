summary.dev.test <-
function(object, ...)
{
	alpha<-0.05
	
	if (!inherits(object, "dev.test")) {
		stop("Error: 'object' must be of class 'dev.test'.")
	}
	
	# Wald Approximation, standard error, confidence interval, coefficient of variation
	#number.of.perm<-function(p.value, cv) { (p.value*(1-p.value))/(cv*p.value)^2 }
	se<-function(p.value, perm=100) { sqrt((p.value*(1-p.value))/perm) }
	ci<-function(p.value, perm=100, alpha=alpha) { 
		tmp<-p.value + cbind(qnorm(1-alpha/2)*se(p.value, perm))%*%c(-1, 1)
		tmp[, 1]<-ifelse(tmp[, 1] < 0, 0 , tmp[, 1])
		tmp[, 2]<-ifelse(tmp[, 2] > 1, 1 , tmp[, 2])
		return(tmp)
	}
	cv<-function(p.value, perm=100) { se(p.value, perm)/p.value }

	x<-unclass(object)
	dat<-as.data.frame(x$test)
	b<-x$number.of.permutations
	
	cat(paste("   full formula: ", paste(as.character(x$formula.full), collapse=" "), "\n", sep=""))
	cat(paste("reduced formula: ", paste(as.character(x$formula.red), collapse=" "), "\n\n", sep=""))
	cat(paste("method: ", as.character(x$method)), "\n", sep="")
	cat(paste("number of variables: ", as.character(x$number.of.variables)), "\n", sep="")
	cat(paste("number of permutations: ", b, "\n\n", sep=""))
	cat("test: \n")
	print(round(dat, ...))
	if(as.character(x$method) == "permutation") {
		a<-cbind(se(p.value=dat[, 3], perm=b), ci(p.value=dat[, 3], perm=b, alpha=alpha), 
		cv(p.value=dat[, 3], perm=b), dat[, 1]/dat[, 2])
		rownames(a)<-rownames(dat)
		colnames(a)<-c("se", "lower.CI", "upper.CI", "coef.of.variation", "deviance/df")
		cat("\n")
		print(round(a, ...))
	}
}
