pcor.test <-
function(x,y,z,use="mat",calcMethod="p",na.rm=TRUE){
	# The partial correlation coefficient between x and y given z
	#
	# pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
	#
	# x and y should be vectors
	#
	# z can be either a vector or a matrix
	#
	# use: There are two methods to calculate the partial correlation coefficient.
	#	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
	#	 Default is "mat".
	#
	# method: There are three ways to calculate the correlation coefficient, 
	#	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
	# 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
	#	    Default is "p".
	#
	# na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
	#        If not, the missing samples will be removed just when the correlation coefficient is calculated.
	#	   However, the number of samples for the p-value is the number of samples after removing 
	#	   all the missing samples from the whole dataset.
	#	   Default is "T".

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z,check.names=FALSE)

	if(use == "mat"){
		p.use <- "Var-Cov matrix"
		pcor = pcor.mat(x,y,z,corMethod=calcMethod,na.rm=na.rm)
	}else if(use == "rec"){
		p.use <- "Recursive formula"
		pcor = pcor.rec(x,y,z,corMethod=calcMethod,na.rm=na.rm)
	}else{
		stop("\'use\' should be either \"rec\" or \"mat\"!\n")
	}

	# print the method
	if(gregexpr("p",calcMethod)[[1]][1] == 1){
		p.method <- "Pearson"
	}else if(gregexpr("s",calcMethod)[[1]][1] == 1){
		p.method <- "Spearman"
	}else if(gregexpr("k",calcMethod)[[1]][1] == 1){
		p.method <- "Kendall"
	}else{
		stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
	}

	# sample number
	n <- dim(na.omit(data.frame(x,y,z,check.names=FALSE)))[1]
	
	# given variables' number
	gn <- dim(z)[2]

	# p-value
	if(p.method == "Kendall"){
		statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  		p.value <- 2*pnorm(-abs(statistic))
	}

	data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use,check.names=FALSE)
}

