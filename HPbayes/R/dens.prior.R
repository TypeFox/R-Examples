dens.prior <-
function(x, pri.lo = c(0, 0, 0, 0, 0, 15, 0, 0), pri.hi = c(0.15, 
    1, 1, 0.25, 15, 55, 0.1, 1.25)) 
{
    y <- (dunif(x[, 1], pri.lo[1], pri.hi[1]) 
    * dunif(x[, 2], pri.lo[2], pri.hi[2]) 
    * dunif(x[, 3], pri.lo[3], pri.hi[3]) 
    * dunif(x[, 4], pri.lo[4], pri.hi[4]) 
    * dunif(x[, 5], pri.lo[5], pri.hi[5]) 
    * dunif(x[, 6], pri.lo[6], pri.hi[6]) 
    * dunif(x[, 7], pri.lo[7], pri.hi[7]) 
    * dunif(x[, 8], pri.lo[8], pri.hi[8])
    	)
    return(y)
}

