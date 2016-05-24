"zemll" <- function(b,x,yo,do) {
n  <- length(yo)
np <- length(b)
r  <- rep(0,n)
f.res <- .Fortran("zemll",
b=to.single(b),x=to.single(x),y=to.single(yo),do=to.integer(do),n=to.integer(n),np=to.integer(np),
r=to.single(r),emll=single(1))
f.res$emll
}


