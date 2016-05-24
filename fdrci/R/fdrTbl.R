fdrTbl <-
function(obs_vec,perm_list,pname,ntests,lowerbound,upperbound,incr=.1,cl=.95,c1=NA){
	# obs_vec is a vector of observed p-values
	# lowerbound and upperbound define -log10(p-value) range over which fdr is computed for a sequence of thresholds
	# If obs_vec and perm_list have high p-values filtered out, then lowerbound should be >= filtering threshold on -log10(p-value) scale 
	# perm_list is a list of dataframes, each of which include a column with permutation p-values
	# pname contains a string that is the name of the permutation p-value column
	# ntests is the number of tests conducted (if no filtering is done, ntests == length(obs_vec)

	plotdat = as.data.frame(matrix(NA,nrow=0,ncol=8))
	names(plotdat) = c("threshold","fdr","ll","ul","pi0","odp","S","Sp")
	thres_ = seq(lowerbound,upperbound,incr)
	for(i in 1:length(thres_)){
	   plotdat[i,"threshold"] = thres_[i]
	   thr = 10^-(thres_[i])
	   tmp = fdr_od(obs_vec,perm_list,pname,ntests,thr,cl=.95,c1)
	   if(length(tmp) == length(plotdat) - 1) plotdat[i,-1] = tmp
	}
	
	return(plotdat)
} # End fdrTbl
