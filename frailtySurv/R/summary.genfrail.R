summary.genfrail <- function(object, ...) {
  dat <- object  
  s <- append(attributes(dat), list(
      n.obs=length(dat$time),
      n.clusters=length(unique(dat$family)),
      mean.cluster=mean(table(dat$family)),
      censor.rate=1 - sum(dat$status)/length(dat$time)
  ))
  
  class(s) <- "summary.genfrail"
  
  s
}