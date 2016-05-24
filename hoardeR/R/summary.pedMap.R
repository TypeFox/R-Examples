summary.pedMap <- function(object, ...){

   pedMapTable <- table(t(as.vector(object$geno)))

   cat("Summary of pedMap object\n")
   cat("---------------\n")
   cat("Map Matrix dimensions:",nrow(object$map),"x",ncol(object$map),"\n")
   cat("Familie Matrix dimensions:",nrow(object$fam),"x",ncol(object$fam),"\n")
   cat("Genotype Matrix dimensions:",nrow(object$geno),"x",ncol(object$geno),"\n")
   cat("No. Allele 'A/A':", pedMapTable[1],"\n")
   cat("No. Allele 'A/B':", pedMapTable[2],"\n")
   cat("No. Allele 'B/B':", pedMapTable[3],"\n")
   cat("No. Allele '0/0':", pedMapTable[4],"\n")
   cat("No. Monomorphic Loci:", object$monomorphic,"\n")
   invisible(object)
} 
