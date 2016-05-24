print.summary.haplin <- function(x, digits = NULL, ...){
#
## SET DIGITS
	if(!missing(digits)) {
		.dig <- .dig.overall <- digits
	}
##	else .dig.overall <- options()$digits
	else {.dig.overall <- 4
		.dig <- 3
	}

.info <- x$info
.alleles <- x$alleles
.nas <- .info$haplos$alleles.nas
.ntri.seq <- x$ntri.seq
#
cat("----Arguments supplied to haplin in this run:----\n")
cat("\n")
print(.info, args.only = T)
cat("\n")
#
cat("----Data summary:----\n")
#
## ACCOUNTING FOR LOST DATA
.summ.tri.tab <- matrix("", ncol = 3, nrow = 4, dimnames = list(rep("", 4), rep("", 3)))
.summ.tri.tab[,1] <- format(c("Cause of loss", "Missing data", "Mendelian incons.", "Unused haplotypes"))
.summ.tri.tab[,2] <- format(c("Triads removed", .ntri.seq[1] - .ntri.seq[2], .ntri.seq[2] - .ntri.seq[3], .ntri.seq[3] - .ntri.seq[4]), justify = "right")
.summ.tri.tab[,3] <- format(c("Triads remaining", .ntri.seq[2], .ntri.seq[3], .ntri.seq[4]), justify = "right")

cat("\nNumber of triads in original file:", .ntri.seq[1], "\n")
cat("\nAccounting for possible loss of triads:\n")
print(.summ.tri.tab, quote = F)
cat("\nTriads remaining for analysis:", .ntri.seq[4], "\n")


cat("\nNOTE: In the following, the most frequent allele\n is printed as upper-case, all others are lower-case\n")

for (i in seq(along = .alleles)){
	cat("\nMarker ", names(.alleles)[i], ":\n", sep = "")
	#
	## PREPARE NUMBERS FOR PRINTING
	.all <- .alleles[[i]]
	.tot <- sum(.all)
	.pros <- round(.all/sum(.all) * 100, 1)
	.all <- c(.all, total = .tot)
	.navn <- names(.all)
	.pros <- c(.pros, total = 100)
	.pr.all <- cbind(.navn, format(.all), format(.pros))
	dimnames(.pr.all) <- list(rep("", dim(.pr.all)[[1]]), c("Allele", "Frequency", "Percent"))
	#
	## PRINT ALLELE FREQUENCIES AND HWE TEST
	#cat("Missing alleles: ", 2*x$HWE.res[[i]]$n.miss.geno, "\n")
	cat("Missing alleles: ", .nas[i], "\n")
	print(.pr.all, quote = F)
	cat("Chi-squared test for HWE, p-value: ", format(x$HWE.res[[i]]$p.value, digits = .dig.overall), "\n")
}# END FOR EACH MARKER
	cat("\n")
#
#
cat("--------\n")
#
.haplo <- x$selected.haplotypes
#
#
cat("Haplotypes removed because of low frequencies:\n")
if (length(.haplo[!.haplo]) == 0) cat("None\n")
else cat(names(.haplo)[!.haplo], "\n")
cat("\n")
#
#
cat("Haplotypes used in the analysis, with coding:\n")
.tmp <- seq(along = .haplo[.haplo])
names(.tmp) <- names(.haplo[.haplo])
print(.tmp, ...)
cat("\n")
#
#
#
## PRINT ESTIMATED PARAMETERS ETC.
cat("----Estimation results:----\n")
if(.info$model$scoretest == "only"){
	cat('NOTE: full estimation results not computed when scoretest == "only"\n\n')
}else{
	print(x$summary.tri.glm, digits = .dig, haplos = names(.haplo)[.haplo], ...)	
}
#
## PRINT LIKELIHOOD RATIO AND/OR SCORE TEST RESULT
cat("\nOverall test for difference between null model (no effects) and full model:\n")
if(.info$model$scoretest %in% c("yes", "no")){
	## PRINT RESULT FROM LIKELIHOOD RATIO TEST
	cat("------------\nLIKELIHOOD RATIO TEST:")
	.ut <- matrix("", ncol = 2, nrow = 4, dimnames = list(rep("", 4), rep("", 2)))
	.ut[,1] <- format(c("Loglike null model:", "Loglike full model:", "df:", "Likelihood ratio p-value:"))
	.ut[,2] <- format(x$loglike, scientific = F, digits = .dig.overall, nsmall = .dig.overall)
	print((.ut), quote = F)
}
if(.info$model$scoretest %in% c("yes", "only")){
	## PRINT RESULT FROM SCORE TEST
	cat("------------\nSCORE TEST:")
	.ut <- matrix("", ncol = 2, nrow = 3, dimnames = list(rep("", 3), rep("", 2)))
	.ut[,1] <- format(c("Score chi-squared value:", "Score df:", "Score p-value:"))
	.ut[,2] <- format(c(x$score[c("chisquared", "df", "pval")]), scientific = F, digits = .dig.overall, nsmall = .dig.overall)
	print((.ut), quote = F)
}
cat("\n(NOTE: The test may be sensitive to rare haplotypes)\n")
#
##
return(invisible(x))
}
