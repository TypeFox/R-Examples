Armitage.MCPerm <-
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
	
    N = case_11+case_12+case_22+control_11+control_12+control_22
    R = case_11+case_12+case_22
    S = N-R
	n0=case_11+control_11
    n1=case_12+control_12
	n2=case_22+control_22
	temp2 = n1+2*n2
	temp3 = n1+4*n2
	low=S*R*(N*temp3-temp2*temp2)
	
	# test for true data
	temp1 = case_12 + 2 * case_22
	upp=N*(N*temp1-R*temp2)*(N*temp1-R*temp2)
	obsTrendValue = upp/low
    obsTrendP = pchisq(obsTrendValue,df=1,lower.tail=FALSE)
	
	# permute and test for true data
	perm_case_11 = rhyper(repeatNum, m = R , n = S, k = n0 )
	perm_control_11=n0 - perm_case_11
	count_case = R - perm_case_11
	count_control=S - perm_control_11
	perm_case_12 = rhyper(repeatNum, m = count_case, n=count_control, k=n1 ) 
    perm_case_22=count_case-perm_case_12
	
	temp1=perm_case_12+2*perm_case_22
	upp=N*(N*temp1-R*temp2)*(N*temp1-R*temp2)
	
	trendValue=upp/low
    trendP=pchisq(trendValue,df=1,lower.tail=FALSE)
	
    maxExp = trendValue[trendValue > obsTrendValue ]
    pValue = length(maxExp)/repeatNum
    list(pValue = pValue, obsStatistic = obsTrendValue , obsP = obsTrendP, 
	    permStatistic = trendValue, permP = trendP)
}
