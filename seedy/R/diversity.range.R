diversity.range <-
function(m.rate, runtime, equi.pop, iterations=10, n.points=100, genomelength=100000, 
         bottle.times=0, bottle.size=1, feedback=1000, makeplot=TRUE, area=TRUE, colline="blue", colarea=rgb(0,0,1,0.4), 
         ref.strain=NULL, init.freq=1, libr=NULL, 
         nuc=NULL, ...) {
  diversity <- NULL
  s.times <- unique(seq(from=1, by=floor(runtime/n.points), to=runtime))
  n.points <- length(s.times)
  
  for (i in 1:iterations) {
    cat("Iteration", i, "\n")
    X <- simulatepopulation(m.rate=m.rate, runtime=runtime, equi.pop=equi.pop, init.freq=init.freq, 
                          bottle.times=bottle.times, bottle.size=bottle.size,libr=libr, 
                          nuc=nuc, genomelength=genomelength, full=TRUE, sample.times=1:runtime,
                          feedback=feedback)
    Y <- plotdiversity(X, s.times, makeplot=FALSE, filter=TRUE)
    diversity <- rbind(diversity, Y)
  }
  if (makeplot) {
    plot(NULL, xlim=c(0,runtime), ylim=c(0,max(diversity)), xlab="Time", ylab="Diversity", ...)
    if (!area) {
      for (i in 1:iterations) {
        lines(s.times, diversity[i,], col=rgb(0,0,1,0.5))
      }
    } else {
      polygon(c(s.times, rev(s.times)), c(apply(diversity,2,quantile, 0.025), rev(apply(diversity,2,quantile, 0.975))), 
              col=colarea, border=NA)
      lines(s.times, apply(diversity, 2, mean))
    }
  }
  return(diversity)
}
