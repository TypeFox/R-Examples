xprod <- function(a,b){
	r <- vector("numeric", 3)
	r[1] <- a[2] * b[3] - a[3] * b[2]
    r[2] <- a[3] * b[1] - a[1] * b[3]
    r[3] <- a[1] * b[2] - a[2] * b[1]
return(r)
}