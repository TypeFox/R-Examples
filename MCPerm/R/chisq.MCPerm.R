chisq.MCPerm <-
function(case_11,case_12,case_22,control_11,control_12,control_22,repeatNum=1000)
{
	# if ( case_11<0 || case_11 != round(case_11) ) {
        # stop("'case_11' must be a positive integer.")
    # }
	# if ( case_12<0 || case_12 != round(case_12) ) {
        # stop("'case_12' must be a positive integer.")
    # }
	# if ( case_22<0 || case_22 != round(case_22) ) {
        # stop("'case_22' must be a positive integer.")
    # }
	
    # if ( control_11<0 || control_11 != round(control_11) ) {
        # stop("'control_11' must be a positive integer.")
    # }
	# if ( control_12<0 || control_12 != round(control_12) ) {
        # stop("'control_12' must be a positive integer.")
    # }
	# if ( control_22<0 || control_22 != round(control_22) ) {
        # stop("'control_22' must be a positive integer.")
    # }
	
	# if ( repeatNum<0 || repeatNum != round(repeatNum) ) {
        # stop("'repeatNum' must be a positive integer.")
    # }
	
    genotypeCount = matrix(c(case_11,case_12,case_22,control_11,control_12,control_22),nrow=2,byrow=TRUE)
	# deal with case and control specifed genotype frequecy is 0
	del=0
    delIndex = c()
    for (n in 1:3) {
        if (all(genotypeCount[, n] == 0)) {
		    del=1
            delIndex = c(delIndex, n)
        }
    }
    if (!is.null(delIndex)) {
        genotypeCount = genotypeCount[, -delIndex, drop = FALSE]
    }
	
	# test for true data 
    chiValue = matrix(0, nrow = 1, ncol = repeatNum )
    chiP = matrix(0, nrow = 1, ncol = repeatNum )
    temp = chisq.test(genotypeCount)
    obsChiValue = temp$statistic
    obsChiP= temp$p.value
	
	# permutate true data and testMethod 
	count_11 = case_11 + control_11
    count_12 = case_12 + control_12
    count_22 = case_22 + control_22
    count_case = case_11 + case_12 + case_22
    count_control = control_11 + control_12 + control_22
	perm_case_11 = rhyper(repeatNum, m = count_case, n = count_control, k = count_11 )
    perm_control_11 = count_11-perm_case_11
	
	count_case = count_case - perm_case_11
	count_control = count_control - perm_control_11
	perm_case_12 = rhyper(repeatNum, m = count_case, n = count_control, k = count_12 ) 
	perm_control_12=count_12-perm_case_12
	 
    perm_case_22=count_case-perm_case_12
	perm_control_22=count_control-perm_control_12
	
    for (i in 1:repeatNum) {
	    randCount=matrix(c(perm_case_11[i],perm_case_12[i],perm_case_22[i],
		    perm_control_11[i],perm_control_12[i],perm_control_22[i]), nrow=2,byrow=TRUE)
	    if(del==1){
            randCount = randCount[, -delIndex, drop = FALSE]
	    }
        temp = chisq.test(randCount)
        chiValue[1, i] = temp$statistic
        chiP[1, i] = temp$p.value
	}
	
	# return result
    maxExp = chiValue[chiValue > obsChiValue]
    pValue = length(maxExp)/repeatNum
    list(pValue = pValue, obsStatistic = obsChiValue, obsP = obsChiP, 
	   permStatistic = chiValue, permP = chiP)
}
