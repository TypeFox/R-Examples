Armitage <-
function(case_11,case_12,case_22,control_11,control_12,control_22)
{
     # if ( control_11<0 || control_11 != round(control_11)) {
         # stop("'control_11' must be a positive integer.")
     # }
	 # if ( control_12<0 || control_12 != round(control_12)) {
         # stop("'control_12' must be a positive integer.")
     # }
	 # if ( control_22<0 || control_22 != round(control_22)) {
         # stop("'control_22' must be a positive integer.")
     # }
	
	 # if ( case_11<0 || case_11 != round(case_11)) {
         # stop("'case_11' must be a positive integer.")
     # }
	 # if ( case_12<0 || case_12 != round(case_12)) {
         # stop("'case_12' must be a positive integer.")
     # }
 	 # if ( case_22<0 || case_22 != round(case_22)) {
         # stop("'case_22' must be a positive integer.")
     # }
	 
	 # Armitage¡¯s trend test for the 2x3 genotype table
     N = case_11+case_12+case_22+control_11+control_12+control_22
     R = case_11+case_12+case_22
     S = N-R
     temp1 = case_12 + 2 * case_22
     n1=case_12+control_12
	 n2=case_22+control_22
	 temp2 = n1+2*n2
	 temp3 = n1+4*n2
	 upp=N*(N*temp1-R*temp2)*(N*temp1-R*temp2)
	 low=S*R*(N*temp3-temp2*temp2)

     statistic = upp/low
     pValue = pchisq(statistic, df = 1, lower.tail = FALSE)
     return(list(statistic = statistic, pValue = pValue))
}
