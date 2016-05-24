print.suest <- function(x, verbose = F, ...){

.suest <- x


cat("\nIndividual score p-values:\n")
#print(.suest$pval.obs)
print(.suest$pval.alt)
cat("\nMultiple testing corrected p-value, based on minimum of individual values:\n")
#print(.suest$pval.obs.corr)
print(.suest$pval.alt.corr)
cat("\nNumber of data lines used in each estimation:\n")
print(.suest$lines.account)

if(verbose){

	cat("\n###################\n")
	cat("Stuff used for testing, ignore:\n")
	print(str(.suest[c("bonferroni", "kill", "pval.alt", "pval.alt.corr")]))
}

.test <- cbind("Individual score p-values:", .suest$pval.alt)
#.test <- cbind("Individual score p-values:", .suest$pval.obs)



dimnames(.test) <- list(rep("", dim(.test)[1]), rep("", dim(.test)[2]))


return(invisible(.test))

}
