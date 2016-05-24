"f.EM.missing" <- function(data, maternal, response = "free", max.EM.iter, x = F, verbose = T, suppressEmWarnings = F, info)
{
## PERFORMS THE ESSENTIALS OF THE EM-ALGORITHM, INCLUDING HANDLING AMBIGUOUS HAPLOTYPES
## AND MISSING GENOTYPE DATA
#
#
## SAKER SOM HAR MED AGGREGERTE DATA AA GJORE
###	.tmpfreq <- rep(1, length(.agg.freq))
###	.agg.freq <- agg.data$freq	#
###	.amb <- data$amb
###	.freqsum <- f.groupsum(X = .freq, INDICES = .amb)
###	.denominator <- f.groupsum(X = rep(1, length(.amb)), INDICES = .amb)	# COULD USE tmp FROM ORIGINAL MATRIX, BUT f.thin HAS BEEN USED IN THE MEANTIME...
# COMPUTE STARTING FREQUENCIES FOR THE EM BY TAKING AVERAGE OVER AMBIGUOUS CATEGORIES:
###	.tmpfreq <- .freqsum/.denominator	#
##	.tmpfreq <- .freq
#
#
#### INITIALIZE EM LOOP: ##################
.tmpfreq <- -1 ## SIGNALS INITIALIZATION
.deviance <- numeric(max.EM.iter)	#
i <- 1	#
.EM.conv <- F
#
## SET UP DESIGN MATRIX, ONE TIME ONLY
.design.matrix <- f.make.design(maternal = maternal, response = response, info = info)
#### EM LOOP: ###########################
repeat{
	# M-STEP:
	#
	## ACTUAL ESTIMATION
	.res <- f.tri.glm(.tmpfreq, design.matrix = .design.matrix, maternal = maternal, info = info)	#
	#
	## PARAMETERS NEEDED TO CHECK CONVERGENCE
	.deviance[[i]] <- .res$result$deviance	# THIS IS TWICE THE MINUS LOG-LIKELIHOOD OF THE GLM
	.coef.new <- .res$result$coefficients
	#
	## PRINT INTERMEDIATE RESULTS (IF REQUESTED):
	if(verbose) {
			cat("EM iter:", sprintf("%-4.0f", round(i)), "|")	#
			cat("GLM deviance:", sprintf("%-12.6g", .deviance[[i]]), "|")
			cat("Coefficients:", sprintf("%-12.6g", .coef.new), "\n")
	}# END if(verbose)
	#
	## STOP WHEN CONVERGED
	if(i > 1) {
		#
		## STOPPING CRITERIA:
		.crit.deviance <- abs(.deviance[i] - .deviance[i - 1]) < 2e-006
		.crit.coef <- max(abs(.coef.new - .coef.old)) < 1e-006
		if(.crit.deviance & .crit.coef){
			.EM.conv <- T
			break
		}
	}# END if(i > 1)
	#
	## BREAK OFF, WITH WARNING, IF max.EM.iter IS EXCEEDED:
	if(i >= max.EM.iter) {
		if(!suppressEmWarnings) warning("Maximum number of EM iterations reached!\n Convergence not yet obtained. Setting max.EM.iter higher may help.")
		break
	}
	#
	## UPDATING VALUES, E-STEP:  
	##
	i <- i + 1
	.coef.old <- .coef.new
	.pred <- .res$pred	#
	.tmpfreq <- f.redistribute(pred = .pred, data = data, info = info)
	#
	## CHECK CONSISTENCY OF NEW VALUES:
		if(any(is.na(.pred))) {
			warning("Missing in predicted values.... There may be a problem with the EM algorithm!")
			.tmpfreq[is.na(.tmpfreq)] <- 1e-005	# FORSIKTIG HER!!!
		}
		if(sum(is.na(.tmpfreq)) > 0) {
			stop("Missing values in EM-updated frequencies!")
		}
next
}# END REPEAT
#
if(x){
	## ADD DESIGN MATRIX IN LAST STEP, IF REQUESTED:
	.res <- f.tri.glm(.tmpfreq, design.matrix = .design.matrix, maternal = maternal, info = info, x = x)#
} # WARNING: REQUIRES .tmpfreq TO BE UNCHANGED AFTER LAST COMPUTATION OF .res!
#
attr(.res, "iter.used") <- i
attr(.res, "EM.conv") <- .EM.conv
#
return(.res)
}
