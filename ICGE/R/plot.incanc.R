plot.incanc <- function(x, ...){
# x of class incanc
   method <- x$method
   
   index <- x$INCAindex                                
   K <- length(index)
   plot(1:K, index , type="b", xlab="Number of clusters", ylab="INCA index",
        sub=paste("Clustering: ", method), xlim=c(1.5, K+0.5), ylim=c(0,1), col="blue")
   cond1 <- !is.null(x$noise[[2]])
   cond2 <- length(x$noise[[2]])!=0
   if (cond1 | cond2){    
     lines(1:K, x$noise[[1]], type="b", col="red")
     legend("topright", legend=c("with noise", "without noise"), pch=21, col=c("red", "blue"))
   }
}
