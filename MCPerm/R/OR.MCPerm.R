OR.MCPerm <-
function(case_allele1, case_allele2, control_allele1, control_allele2,repeatNum=1000)
{
    # if (case_allele1<0 || case_allele1 != round(case_allele1)) {
         # stop("'case_allele1' must be a positive integer.")
    # }
	# if (case_allele2<0 || case_allele2 != round(case_allele2)) {
         # stop("'case_allele2' must be a positive integer.")
    # }
	# if (control_allele1<0 || control_allele1 != round(control_allele1)) {
         # stop("'control_allele1' must be a positive integer.")
    # }
	# if (control_allele2<0 || control_allele2 != round(control_allele2)) {
         # stop("'control_allele2' must be a positive integer.")
    # }
	
	# if ( repeatNum<0 || repeatNum != round(repeatNum) ) {
        # stop("'repeatNum' must be a positive integer.")
    # }	
    
	# test for true data 
	OR1=(case_allele1*control_allele2)/(case_allele2*control_allele1)
    OR2=1/OR1
    # Determine the risk allele(OR>1)
	if(OR1>OR2){
	     risk_allele=1
		 OR=OR1
	}else{
	     risk_allele=2
		 OR=OR2
	}
	
    allele1=case_allele1+control_allele1
	allele2=case_allele2+control_allele2
	case_count=case_allele1+case_allele2
	control_count=control_allele1+control_allele2
	
	perm_case_allele1=rhyper(repeatNum,case_count,control_count,allele1)
    perm_control_allele1=allele1-perm_case_allele1
	perm_case_allele2=case_count-perm_case_allele1
    perm_control_allele2=control_count-perm_control_allele1
	if(risk_allele==1){
		 ORValue=(perm_case_allele1*perm_control_allele2)/(perm_case_allele2*perm_control_allele1)
	}else{
		 ORValue=(perm_case_allele2*perm_control_allele1)/(perm_case_allele1*perm_control_allele2)
	}
	
    maxExp = ORValue[ORValue > OR]
    pValue = length(maxExp)/repeatNum
    list(risk_allele=risk_allele,pValue = pValue, obsOR = OR,permOR = ORValue)
}
