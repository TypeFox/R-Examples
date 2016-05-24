"FunnelPlot" <-
function(x)
 {
 rxy <- x$Rxy
 N <- x$n
 rb <- rbar(x)
 plot(rxy,N, xlab="Effect Sizes", ylab="Sample Sizes", main="Funnel Plot")
 abline(v=rb)
 }

