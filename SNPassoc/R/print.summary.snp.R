print.summary.snp<-function(x,...)
{
       cat("Genotypes: \n")
       print(x$genotype.freq, na="",...)
       cat("\n")
       cat("Alleles: \n")
       print(x$allele.freq, na="", ...)
       cat("\n")
       cat("HWE (p value):", x$HWE, "\n")
}
