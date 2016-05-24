f.preliminary.freq.new <- function(data, info){
##
## COMPUTES PRELIMINARY FREQUENCY DISTRIBUTION OF HAPLOTYPES, 
## BY A SIMPLE EM COMPUTATION
## 
## OUTPUTS A FREQUENCY VECTOR FITTING IN SIZE FOR data (PREDICTED FREQUENCIES),
## WITH ATTRIBUTES CONTAINING A HAPLOTYPE FREQUENCY TABLE
##
## IMPORTANT: ASSUMES m1, m2, f1, f2 (OR c1, c2) TO HAVE NUMERIC CODING! 
## AND THAT HAPLOTYPES ARE ORDERED IN THE "NATURAL" ORDER, WITH THE FIRST 
## LOCUS COUNTING MOST QUICKLY,
## ETC, E.G. 1-1, 2-1, 3-1, 1-2, 2-2, ....
##
## NOTE: DOES *NOT* SET NAMES FOR THE HAPLOTYPES IN THE OUTPUT
##
#
## PREPARATIONS:
#
.ntri <- length(unique(data$ind))
design <- info$model$design
max.EM.iter <- info$control$max.EM.iter
verbose <- info$control$verbose
.xchrom <- info$model$xchrom
#
## CREATE A .freq VECTOR WITH UNIFORM DISTRIBUTION WITHIN FAMILIES:
	.freq <- 1/f.groupsum(X = rep(1, length(data$ind)), INDICES = data$ind)
#
#
if((design == "triad" | design == "cc.triad") & !.xchrom){
## LIST ALL HAPLOTYPES IN SINGLE VECTOR:
	.all <- c(data$m1, data$m2, data$f1, data$f2)
	.ncols <- 4
#
}
if((design == "triad" | design == "cc.triad") & .xchrom){
## LIST ALL HAPLOTYPES IN SINGLE VECTOR:
	.all <- c(data$m1, data$m2, data$f1)
	.ncols <- 3
#
}
#
if(design == "cc"){
## LIST ALL HAPLOTYPES IN SINGLE VECTOR:
	.all <- c(data$c1, data$c2)
	.ncols <- 2
#
}
.all.freq <- rep(.freq, .ncols)
#
## JUST FOR SAFETY'S SAKE:
if(is.factor(.all) | is.character(.all)) stop("Problem with data type in f.preliminary.freq...")
### if((design != "triad") & (design != "cc")) stop()

#
## PREPARE FOR EM ALGORITHM:
#
if(verbose) cat("\nRunning EM for preliminary estimates of haplotype frequencies...  ")
#
i <- 0
.oldtable <- 0
.epsilon <- 1e-5 # SHOULD BE SUFFICIENT FOR MAX OF CHANGE IN REL. FREQUENCIES
.alleles <- info$haplos$alleles
.haplotypes <- 1:(prod(sapply(.alleles, length)))
.aux.haplotypes <- seq(along = .haplotypes)
.aux.freq <- rep(0, length.out = length(.haplotypes))
#
#
## LOOP OVER EM:
	repeat{
		i <- i + 1
		if(i >= max.EM.iter) {
			warning("Maximum number of EM iterations reached!\n Convergence not yet obtained. Setting max.EM.iter higher may help.")
			break}
#
#		FREQUENCY TABLE OVER ALL HAPLOTYPES:		
		.table <- tapply(c(.aux.freq, .all.freq), c(.aux.haplotypes, .all), sum) ## THE .aux ENSURES ALL HAPLOTYPES ARE PRESENT
		.table <- .table/sum(.table)
		if(verbose & F) print(round(.table, 8))
#
#		STOP WHEN FREQUENCY TABLE HAS CONVERGED SUFFICIENTLY:
		if(max(abs(.table - .oldtable)) < .epsilon) break
		.oldtable <- .table
#
#		PREDICT FROM MULTIPLICATIVE MODEL:
		.pred <- matrix(.table[.all], ncol = .ncols)
		#
		if((design == "triad" | design == "cc.triad") & !.xchrom){
			.pred <- .pred[,1] * .pred[,2] * .pred[,3] * .pred[,4] * .ntri # AVOID apply FOR SPEED WITH LARGER MATRICES
		}
		if((design == "triad" | design == "cc.triad") & .xchrom){
			.pred <- .pred[,1] * .pred[,2] * .pred[,3] * .ntri # AVOID apply FOR SPEED WITH LARGER MATRICES
		}
		if(design == "cc"){
#			if(.xchrom) stop ("Not implemented!")
			.pred <- .pred[,1] * .pred[,2] * .ntri # AVOID apply FOR SPEED WITH LARGER MATRICES
		}
#
#		RENORMALIZE WITHIN EACH FAMILY:
		.predsum <- f.groupsum(X = .pred, INDICES = data$ind)
		.tmpfreq <- ifelse(.predsum > 0, .pred/.predsum, 0)
#
		.all.freq <- rep(.tmpfreq, .ncols)
} # END repeat
#
#
.freq.ut <- .tmpfreq
attr(.freq.ut, "prelim.haplotype.freq") <- .table
#
if(abs(sum(.freq.ut) - .ntri) > 0.001 * .ntri) stop("Potential problem in preliminary EM frequencies!")	# JUST CHECKING THAT "MASS" DOES NOT DISAPPEAR
#
if(verbose) cat("Done\n")
#
return(invisible(.freq.ut))
}
