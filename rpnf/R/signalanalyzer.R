#' analyze transitions of signal states
#'
#' @param signal Signal to identify
#' @param probability Report probability
signalanalyzer <- function(signal,probability=TRUE) {
  t2<-table(signal[1:(length(signal)-1)],signal[2:length(signal)])
  t2m<-as.matrix(t2)
  for (i in 1:nrow(t2m)) t2m[i,i]=0
  if (probability)
    for (i in 1:nrow(t2m)) t2m[i,]=t2m[i,]/sum(t2m[i,])
  # replace 0 with NA
  for (i in 1:length(t2m))
    if (t2m[i]==0)
      t2m[i]<-NA
  t2m
}
