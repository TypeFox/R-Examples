permuteGenotypeCount <-
function (case_11, case_12, case_22, control_11, control_12, control_22,n) 
{
     # if( case_11<0 || case_11 != round(case_11) ){
	     # stop("case_11 should be a positive integer.")
	 # }
	 # if( case_12<0 || case_12 != round(case_12) ){
	     # stop("case_12 should be a positive integer.")
	 # }
	 # if( case_22<0 || case_22 != round(case_22) ){
	     # stop("case_22 should be a positive integer.")
	 # }
	 
	 # if( control_11<0 || control_11 != round(control_11) ){
	     # stop("control_11 should be a positive integer.")
	 # }
	 # if( control_12<0 || control_12 != round(control_12) ){
	     # stop("control_12 should be a positive integer.")
	 # }
	 # if( control_22<0 || control_22 != round(control_22) ){
	     # stop("control_22 should be a positive integer.")
	 # }
	 
     count_11 = case_11 + control_11
     count_12 = case_12 + control_12
     count_22 = case_22 + control_22
     count_case = case_11 + case_12 + case_22
     count_control = control_11 + control_12 + control_22
	 
     perm_case_11 = rhyper(n, m = count_case , n = count_control, k = count_11 )
     perm_control_11 = count_11-perm_case_11
	 
	 count_case = count_case - perm_case_11
	 count_control = count_control - perm_control_11
	 perm_case_12 = rhyper(n, m = count_case, n=count_control, k=count_12 ) 
	 perm_control_12=count_12-perm_case_12
	 
     perm_case_22=count_case-perm_case_12
	 perm_control_22=count_control-perm_control_12
     
     result=list('perm_case_11'=perm_case_11,'perm_case_12'=perm_case_12,'perm_case_22'=perm_case_22,
	   'perm_control_11'=perm_control_11,'perm_control_12'=perm_control_12,'perm_control_22'=perm_control_22)
     return(result)
}
