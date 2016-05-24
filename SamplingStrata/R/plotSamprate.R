plotSamprate <- function (solution, dom) {
  df <- as.data.frame(solution$aggr_strata[solution$aggr_strata$DOM1 == dom,])
  plot(df$SOLUZ/df$N,
       ylab="sampling units / total units",
       xlab="strata",
	   type="h")
  title(paste("Sampling rate per stratum in domain ",dom,sep=""))
}    