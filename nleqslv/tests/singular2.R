# http://stackoverflow.com/questions/29134996/solving-nonlinear-equation-in-r

# wants to know if system has closed form solution
# I want to see how nleqslv behaves

set.seed(29)

library(nleqslv)

print.result <- function(z, do.print.xf=FALSE) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

f <- function(X, a, b, c1, c2, c3) {
    Y <- numeric(3)
    x <- X[1]
    y <- X[2]
    z <- X[3]
    Y[1] <- x + y - x*y - c1
    Y[2] <- x + z - x*z - c2
    Y[3] <- a*y + b*z - c3
    return(Y)
}

Jac <- function(X, a, b, c1, c2, c3) {
    J <- matrix(0,nrow=3,ncol=3)
    x <- X[1]
    y <- X[2]
    z <- X[3]

    J[1,1] <- 1-y
    J[2,1] <- 1-z
    J[3,1] <- 0
    J[1,2] <- 1-x
    J[2,2] <- 0
    J[3,2] <- a
    J[1,3] <- 0
    J[2,3] <- 1-x
    J[3,3] <- b
    J
}

a <- 1
b <- 1
c1 <- 2
c2 <- 3
c3 <- 4

# exact solution
x <- (a*c1+b*c2-c3)/(a+b-c3)
y <- (b*c1-b*c2-c1*c3+c3)/(-a*c1+a-b*c2+b)
z <- (a*(c1-c2)+(c2-1)*c3)/(a*(c1-1)+b*(c2-1))
xsol <- c(x,y,z)

X.start <- c(1,2,3)
z1 <- nleqslv(X.start,f,Jac,a=a,b=b,c1=c1,c2=c2,c3=c3,
                method="Newton",control=list(trace=0,allowSingular=TRUE))

z2 <- nleqslv(X.start,f,Jac,a=a,b=b,c1=c1,c2=c2,c3=c3,
        method="Broyden",control=list(trace=0,allowSingular=TRUE))

all.equal(z1$x,xsol)
all.equal(z2$x,xsol)
print.result(z1)
print.result(z2)
