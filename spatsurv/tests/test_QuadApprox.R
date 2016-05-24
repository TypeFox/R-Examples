library(spatsurv)

func <- function(par){
    x <- par[1]
    y <- par[2]
    z <- par[3]
    
    return(0 + x + 2*y + 3*z + 2*x^2 + y^2 + 2*z^2 + 3*x*y + 4*x*z + 5*y*z + rnorm(1))
}

ans <- QuadApprox(fun=func,npts=10,argRanges=list(c(-5,5),c(-5,5),c(-5,5)),plot=FALSE)
print(ans)

# ans1 <- quadapprox(fun=func,xseq=seq(-5,5,length.out=20),yseq=seq(-5,5,length.out=20))
# print(ans1)