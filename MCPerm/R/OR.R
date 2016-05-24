OR <-
function (case_allele1, case_allele2, control_allele1, control_allele2) 
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
	
	# odd ratio
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
	
	# 95%CI (Confidence Interval)
    LnOR = log(OR)
    SELnOR = sqrt(1/case_allele1 + 1/case_allele2 + 1/control_allele1 + 1/control_allele2)
    lowLnOR = LnOR - 1.96 * SELnOR
	upLnOR = LnOR + 1.96 * SELnOR
    lowerOR = exp(lowLnOR)
    upperOR = exp(upLnOR)
    CI = matrix(c(lowerOR, upperOR), nrow = 1)
	
	# return result
    return(list(risk_allele=risk_allele,OR = OR, CI = CI))
}
