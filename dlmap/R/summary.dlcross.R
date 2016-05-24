`summary.dlcross` <- 
function(object, ...)
{
  cat("\nThis is an object of class dlcross.\nSummary of genetic and phenotypic data:\n\n")

  gen <- as.matrix(object$dfMerged[, -c(1:object$nphe)])
  if ((type <- attr(object, "type"))=="other")
  cat("\t This is an association mapping population.\n\n") else {
  cat("\t This is a", type, "population.\n")}
  cat("\t No. individuals:         ", ngen(object),"\n")
  cat("\t No. phenotypic traits:   ", object$nphe,"\n")
  cat("\t Percent phenotyped:      ", apply(object$dfMerged[,1:object$nphe], 2, function(x) return(sum(!is.na(x))/length(x)*100)),"\n", fill=TRUE)
  cat("\t No. chromosomes:         ", length(object$map), "\n")
  cat("\t Total markers:           ", sum(nmrk(object)),"\n")
  cat("\t Percent genotyped:       ", 100*sum(!is.na(gen))/length(as.vector(gen)),"\n")
  cat("\t No. markers per chr:     ", nmrk(object),"\n", fill=TRUE)
  cat("There are", ngen(object),"unique genotypes and", nphen(object), "unique phenotyped individuals.\n", fill=TRUE)
}
