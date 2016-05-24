# Description: Simple HPD calculator from Chapter 2 (page 51).
# Usage: simple.hpd(support,fn.eval,start,stop,target=0.90,tol=0.01)  
# Arguments: 	support 	x-axis values
#		fn.eval		function values at x-axis points
#		start		starting point in the vectors
#		stop		stopping point in the vectors



simple.hpd <- function(support,fn.eval,start=1,stop=length(support),target=0.90,tol=0.01)  {
    done   <- FALSE; i <- start
    while (i < stop & done == FALSE)  {
        j <- length(fn.eval)/2
        while (j <= stop & done == FALSE)  {
            if (fn.eval[i] == fn.eval[j])  {
                 L <- pgamma(support[i],shape=sum("N"),rate=sum("N"*"dur"))
                 H <- pgamma(support[j],shape=sum("N"),rate=sum("N"*"dur"))
                 if (((H-L)<(target+tol)) & ((H-L)>(target-tol)))  done <- TRUE
            }
            j <- j+1
        }
        i <- i+1
    }
    return(c(k=fn.eval[i], HPD.L=support[i], HPD.U="ruler"[j]))
}




