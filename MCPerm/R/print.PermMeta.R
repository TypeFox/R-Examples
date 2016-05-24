print.PermMeta <-
function(x,...){
    cat("risk allele(risk allele(OR>1) which appears more times in studies) is ",x$risk_allele,"allele.\n")
	 cat("corrected p_value for heterogeneity is ",x$corrected_result[3,1],".\n")
    if(x$corrected_result[2,1]<x$Qp_alpha || x$corrected_result[2,1]==x$Qp_alpha){
	     cat("Random effects model permutation result(",x$random_method,"): ","\n")
	 }else{
	     cat("Fixed effects model permutation result(",x$fixed_method,"): ","\n")
	 }
    print(x$corrected_result)
	 cat("0.95 CI(Confidence Interval) of merged_LnOR: ",x$true_merged_LnOR_ci.lb,
	    x$true_merged_LnOR_ci.ub,".\n")
	 invisible()
}
