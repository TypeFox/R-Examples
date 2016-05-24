f.final.loglike <- function(data, pred, info, type = "EM"){
##
## COMPUTES THE MAXIMUM LOG-LIKELIHOOD (UP TO A CONSTANT) FOR THE FINAL RESULT
## NOTE: TAKES INTO ACCOUNT MISSING INFORMATION, I.E. IS not THE FULL
## LIKELIHOOD IN EM BUT RATHER THE CORRECT MAXIMUM LIKELIHOOD FROM THE OBSERVED
## DATA
##
## data IS A DATA FRAME WITH A TRIAD INDICATOR SHOWING ALL POSSIBLE HAPLOTYPE
## COMBINATIONS FOR THAT TRIAD, pred ARE THE PREDICTED FREQUENCIES IN THE 
## MAXIMIZED FULL LIKELIHOOD, FOR ALL HAPLOTYPE COMBINATIONS IN STANDARD
## ORDERING
##
## type = "EM" IS THE ONE TO USE. type = "full" IS (UP TO A CONSTANT) THE 
## SAME AS THE ONE OBTAINED FROM -res$result$deviance/2, I.E. THE GLM 
## LOGLIKE. 
## SO THAT ONE GETS THE SAME RESULT FROM 
## -res.0$result$deviance + res$result$deviance
## AS FROM
## 2*(.loglike.0 - .loglike)
## WHEN
## .loglike.0 <- f.final.loglike(data = data, pred = res.0$pred, info = info, type = "full")
## .loglike <- f.final.loglike(data = data, pred = res$pred, info = info, type = "full")
#
## STANDARDIZE TO PROBABILITIES, AS IN A MULTINOMIAL:
.prob <- pred/sum(pred) ## (NOT REALLY NEEDED SINCE ONLY A CONSTANT?)
#
##
## MATCH PREDICTED PROBABILITIES TO ORIGINAL DATA:
.pos <- f.pos.match(data = data, info = info)
.prob <- .prob[.pos]
#
##
if(type == "full"){
	## CORRESPONDS TO THE FULL, UNCORRECTED GLM LOGLIKE, UP TO A CONSTANT
	.probsum <- f.groupsum(.prob, data$ind)
	.probnorm <- .prob/.probsum
	.loglike <- sum(.probnorm * log(.prob))
}
if(type == "EM"){
	## THIS IS THE ONE TO USE, CORRECTS FOR EM UNCERTAINTY
	#
	## SUM PREDICTED PROBABILITIES OVER AMBIGUITIES FOR EACH TRIAD:
	.prob <- tapply(.prob, data$ind, sum)
	#
	## COMPUTE LOG-LIKELIHOOD:
	.loglike <- sum(log(.prob)) # NOTE: FREQUENCIES ARE 1
}
#	
## FINISH:
return(.loglike)
}
