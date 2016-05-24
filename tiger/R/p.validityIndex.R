p.validityIndex <- function(result, validity.max){
 if(!is.na(validity.max)){
        xmax <- 1
        while(any(result$validity[xmax:length(result$validity)]<validity.max) & xmax<result$maxc) xmax <- xmax+1
        if(xmax!=result$maxc){
		ylim=c(0,validity.max)
        }
    } else {
        xmax=length(result$validity) + 1
    }

   if(!exists("ylim")) ylim=c(0,max(result$validity, na.rm=TRUE))

   plot(result$validity[1:xmax], ylab=expression(V[XB]), xlab="Number of clusters", type="b", lty=2, ylim=ylim)
}
