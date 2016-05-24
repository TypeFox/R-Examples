plotsnpfreq <- function(data, timepoint=1, type="S", ...) {

  obs.freq <- data$obs.freq[[timepoint]]
  obs.strain <- data$obs.strain[[timepoint]]
  librstrains <- data$librstrains
  
  mutfreq <- numeric(max(unlist(data$libr)))
  
  for (i in 1:length(obs.strain)) {
    marker <- which(librstrains==obs.strain[i])
    if (marker>1) {
      mutfreq[data$libr[[marker]]] <- mutfreq[data$libr[[marker]]]+obs.freq[i]
    }
  }
  
  no.sites <- NULL
  x <- sort(unique(mutfreq))
  x <- x[-which(x==0)]
  if (length(x)>0) {
    for (i in 1:length(x)) {
      no.sites <- c(no.sites, sum(mutfreq>=x[i]))
    }
    x <- c(x, x[length(x)])
    no.sites <- c(no.sites, 0)
    mutant.frequency <- x/sum(obs.freq)
    plot(mutant.frequency, no.sites, xlim=c(0,1), type=type, ...)
  } else {
    warning("No polymorphic sites!")
    plot(NULL, xlim=c(0,1), ylim=c(0,100), ...)    
  }
  return(invisible(rbind(mutant.frequency,no.sites)))
}