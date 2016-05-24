##  WARNINGS: All the codes below are for Exiqon LNA microRNA
##  microarray data analysis. Don't generalize these codes for other
##  microarray platforms.


print.exiqon <- function(x, ...){
  cat("\nClass type: 'exiqon' ( ", nrow(x[[1]]), " obs.)")
  cat("\n")
  print(summary(as.data.frame(x[[1]][c("Signal","Background")])),...)
  invisible(x)
}

summary.exiqon <- function(object, ...){
  cat("\nClass type: 'exiqon'")
  cat("\nSamples:\n")
  print(levels(as.factor(object[[1]]$Sample)))
  cat("\nTreatments:\n")
  print(levels(as.factor(object[[1]]$Treat)))
  genes = levels(object[[1]]$Gene)
  sele = tolower(substr(genes,1,3)) == 'hsa'
  genes = genes[sele]
  cat("\nTotal hsa-miRs: ", length(genes), "\n\n")
  invisible(object)
}

