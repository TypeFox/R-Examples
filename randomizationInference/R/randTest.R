#conduct randomization test and calculate randomization-based p-value
#input: outcomes (y), assignments (w), desired number of permutations (nrand),
#	function to calculate test statistic (calcTestStat),
#	function to calculate potential outcomes under the null hypothesis (calcPO)
#	calcOptions=list of options for calculating potential outcomes and/or test statistic (if necessary),
#	randOptions=list of type/block/row/col options for randomization (if necessary, defaults to complete randomization),
#	userRand=function for user-defined randomization scheme
#	userOptions=list of options for user-defined randomization scheme
#	effect under the null hypothesis (null)
#	alternative hypothesis (defaults to "two.sided")
#output: vector of test statistics based on permuted assignments, observed test statistic,
#	selected alternative hypothesis, randomization-based p-value

randTest=function(y,w,nrand,calcTestStat,calcOptions=NULL,calcPO=zeroEffect,poOptions=NULL,
			randOptions=list(type=c("complete","block","Latin","user.defined"),
		      block = NULL, row = NULL, col = NULL), userRand = NULL, userOptions = NULL,
			null = NULL, alternative=c("two.sided", "greater", "less")){
	alternative=match.arg(alternative)
	randOptions$type=match.arg(randOptions$type,c("complete","block","Latin","user.defined"))

	#convert w for randomization scheme
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)

	#observed test statistic
	if(is.null(calcOptions)) obs_stat=calcTestStat(y,w)
	else obs_stat=calcTestStat(y,w,calcOptions)

	#each list element is a hypothetical random assignment
	if(randOptions$type=="complete"){
		wperm=completeRand(w,nrand)
	}else if(randOptions$type=="block"){
		wperm=blockRand(w,nrand,randOptions$block)
	}else if(randOptions$type=="Latin"){
		wperm=latinRand(w,nrand,randOptions$row,randOptions$col)
	}else if(randOptions$type=="user.defined"){
		#userRand is a function that takes nrand and userOptions as input
		#and outputs a list of random assignments
		wperm=userRand(nrand,userOptions)
	}

	#potential outcomes and corresponding test statistics for each permutation
	if(is.null(calcOptions)){
		if(is.null(poOptions)){
			perm_stats=sapply(1:nrand, function(j) calcTestStat(calcPO(y,w,wperm[[j]]),wperm[[j]]))
		}else{
			perm_stats=sapply(1:nrand, function(j) calcTestStat(calcPO(y,w,wperm[[j]],poOptions),wperm[[j]]))
		}
	}else{
		if(is.null(poOptions)){
			perm_stats=sapply(1:nrand, function(j) calcTestStat(calcPO(y,w,wperm[[j]]),wperm[[j]],calcOptions))
		}else{
			perm_stats=sapply(1:nrand, function(j) calcTestStat(calcPO(y,w,wperm[[j]],poOptions),wperm[[j]],calcOptions))
		}
	}

	#Randomization-based p-value
	if(is.null(null)) null=rep(0,length(obs_stat))
	if(length(obs_stat)==1){
		if(alternative=="two.sided"){
			pvalue=sum(abs(perm_stats-null)-abs(obs_stat-null)>=0)/nrand
		}else if(alternative=="less"){
			pvalue=sum(perm_stats-obs_stat<=0)/nrand
		}else if(alternative=="greater"){
			pvalue=sum(perm_stats-obs_stat>=0)/nrand
		}
	}else{
		if(alternative=="two.sided"){
			pvalue=rowSums(abs(perm_stats-null)-abs(obs_stat-null)>=0)/nrand
		}else if(alternative=="less"){
			pvalue=rowSums(perm_stats-obs_stat<=0)/nrand
		}else if(alternative=="greater"){
			pvalue=rowSums(perm_stats-obs_stat>=0)/nrand
		}
	}

	#output
	if(is.vector(perm_stats)==FALSE) perm_stats=t(perm_stats)
	list(perm_stats=perm_stats,obs_stat=obs_stat,null=null,alternative=alternative,pvalue=pvalue)
}
