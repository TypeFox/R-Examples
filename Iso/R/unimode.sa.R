unimode.sa <- function(y,lmode) {
# Note that the corresponding "x" vector is taken to be 1:n, so
# "lmode" makes most sense if it is one of 1, 1.5, 2, 2.5, ... n-1,
# n-0.5, n.  It can take other values but.  Results are based on
# size comparisons of y with lmode.
#
n  <- length(y)
x  <- 1:n
y1 <- y[x<lmode]
y2 <- y[x>lmode]
n1 <- length(y1)
n2 <- length(y2)
if(n1 <=1 ) return(pava(y,decreasing=TRUE))
if(n2 <=1 ) return(pava(y))
yh1 <- if(n1>0) pava(y1) else NULL
yh2 <- if(n2>0) pava(y2,decreasing=TRUE) else NULL
if(n1+n2==n) {
    yh  <- c(yh1,yh2)
} else {
    yh2 <- rev(yh2)
    o   <- order(c(yh1,yh2))
    r   <- rank(c(yh1,yh2))
    ys  <- c(c(yh1,yh2)[o],y[n1+1])
    yhs <- pava(ys)
    yyy <- (yhs[-n])[r]
    s1  <- seq(to=n1,length=n1)
    s2  <- seq(to=n-1,length=n2)
    yh <- c(yyy[s1], yhs[n], rev(yyy[s2]))
}
yh
}
