
#Sept 2013: per anando mahto request: return NONinvisible

seqle <- function(x, incr=1L, prec=.Machine$double.eps ^ 0.5){ 

    if(!is.numeric(x)) x <- as.numeric(x) 
    n <- length(x)  
    y <- abs(x[-1L] - x[-n] - incr) > prec
    ii <- c(which(y|is.na(y)),n) 
    foo<- list( lengths = diff(c(0L,ii)),  values = x[head(c(0L,ii)+1L,-1L)])
 # The only method for  class 'rle' is "print", so use it here as well.
   class(foo)<-'rle'
   return(foo)
} 

inverse.seqle<-function(x,incr=1L) {
#error checker stolen from inverse.rle:
if (is.null(le <- x$lengths) || is.null(v <- x$values) ||  length(le) != length(v)) {
	stop("invalid 'seqle' structure")
}
theseq<-vector()
# be aware that, for floats, this may not be an exact reconstruction.
for(jj in 1: length(le) ) {
	theseq <- c(theseq, seq(from=x$values[jj],by=incr,length=x$lengths[jj]) )
}
return(theseq)
}
