critical.ci.values <- function(loci2,locus.empirical2){
# Function used within the Bootstrapping.CI.r function, in order to calculate lower and upper confidence interval boundaries
lower.difference <- abs(mean(as.numeric(as.vector(loci2[,1])),na.rm=TRUE)-quantile(as.numeric(as.vector(loci2[,1])),.025,na.rm=TRUE))
upper.difference <- abs(mean(as.numeric(as.vector(loci2[,1])),na.rm=TRUE)-quantile(as.numeric(as.vector(loci2[,1])),.975,na.rm=TRUE))


critical.value.loci <- c(as.numeric(as.vector(locus.empirical2[,1]))-lower.difference,
as.numeric(as.vector(locus.empirical2[,1]))+upper.difference)
invisible(critical.value.loci)

}
