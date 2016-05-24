#################################################################################################
# this function perform the prune process based on Fishers exact test for the splits of the levels
# that we designate, starting from the bottom level.
# It will return a vector of flag indicating the pruned splits, with flag 1 and not pruned splits
# with flag 0.
############################################################################################################

# Ti              a vector of the observed time Ti of each patient, Ti is either the the observed
#                 failture time, or observed censored time.
# delta           a vector of censoring indicator of each patient, censored=0,uncensored=1
# M               total level of the MRH model used
# NumPrune        The number of levels is subject to pruning, starting from the bottom level. For example, M=5 and NumPrune=3
#                 means the 3, 4 and 5 level of  this five-level MRH model are subject to pruning.
# censortime	    the censoring time in the study
# alpha           type one error in Fisher's exact test

Prune=function(Ti, delta, M, NumPrune, maxStudyTime = NULL, alpha = 0.05){ 
	
	censortime = maxStudyTime

	if(missing(NumPrune)){	NumPrune = M	}
	if(NumPrune <= 0){	stop("The number of pruned levels (NumPrune) must be greater than 0.")	}
	if(length(Ti)!=length(delta)){	stop("A failure time and a censoring indicator are needed for each subject") }
		
	indc=rep(0,2^M-1) 
	Time = Ti
	Censor = delta
	if(is.null(censortime)){	censortime = max(Ti)	}
	binwidth = censortime/2^M
    Censor[Time > censortime] = 0
		
	bin = getbin.and.ratio(censortime, Time[delta == 1], M)$bin+1
	bin[bin > 2^M] = 2^M
	bincount = table(factor(bin, levels = 1:2^M))[]
	bincensor = getbin.and.ratio(censortime, Ti[delta == 0], M)$bin+1
	bincensor[bincensor > 2^M] = 2^M
	censor_bincount = table(factor(bincensor, levels = 1:2^M))[]
	censor_bincount[2^M] = censor_bincount[2^M]-length(which(Time[delta == 0] > censortime))
	censor_bincount = c(censor_bincount, length(which(Time[delta == 0] > censortime)))
		
	# Save observed failures in each bin from level 1 to level M in a vector
	vec_bincount = NULL
	for(i in M:1){ 
		for(j in 1:(2^M/2^(i-1))){	vec_bincount = c(vec_bincount, sum(bincount[1:(2^(i-1))+(j-1)*2^(i-1)]))	}
	}
	
	# save observed censored patients in each bin from level 1 to level M in a vector, not include the
	# counts get censored when the study ends
	vec_censor_bincount = NULL
	for(i in M:1){
		for(j in 1:(2^M/2^(i-1))){	vec_censor_bincount = c(vec_censor_bincount, sum(censor_bincount[1:(2^(i-1))+(j-1)*2^(i-1)]))	}
	}
	
	# save patients at risk at the beginning in each bin from level 1 to level M in a vector
    vec_risk = NULL
	totalFails = vec_censor_bincount + vec_bincount
	for(i in 1:M){
		vec_risk = c(vec_risk, length(Time)-c(0, cumsum(totalFails[1:(2^i)+sum(2^(0:(i-1)))-1]))[-(2^i+1)])
	}
										
	# use Fishers exact test to decide the splits which are needed to be pruned and flag them with 1; otherwise
	# flag them with 0.
    flag=NULL
    for( i in 1: (length(vec_bincount)/2)){
		if(vec_bincount[2*i-1]==0  &&  vec_bincount[2*i]==0){	flag=c(flag,1)
		} else if(vec_bincount[2*i-1]==0 || vec_bincount[2*i]==0){	flag=c(flag,0)
		} else {
			XX=matrix(c(vec_bincount[2*i-1], vec_bincount[2*i],
						vec_risk[2*i-1]-vec_bincount[2*i-1], vec_risk[2*i]-vec_bincount[2*i]),2,2,byrow=TRUE)
			temp=fisher.test(XX)
			flag=c(flag,temp[[1]])}
    }

	# Make the pruning indicator, which equals 1 (no pruning) if the p-value from Fishers exact test is > alpha
    indc[which(flag>alpha)]=1

	# If the user only wants a certain number of bins pruned, set all lower bins to 0.
    prune_indc = rep(0,2^M-1) 
    for(mval in M:(M-NumPrune+1)){
		for(pval in 0:(2^(mval-1)-1)){	prune_indc[2^(mval-1)+pval] = indc[2^(mval-1)+pval]	}
	}
    
    return(prune_indc)
	
}


