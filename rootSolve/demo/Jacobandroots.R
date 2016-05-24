

example(uniroot.all)
example(gradient)

# Demonstration of how the banded and full jacobian look like...

mod <- function (t=0,y=c(1,2,3,4), parms=NULL,...)
{
 dy1<-y[1] + 2*y[2]
 dy2<-3*y[1] + 4*y[2] + 5*y[3]
 dy3<-6*y[2] + 7*y[3] + 8*y[4]
 dy4<-9*y[3] + 10*y[4]
 return(as.list(c(dy1,dy2,dy3,dy4)))
}

jacobian.band(y=c(1,2,3,4),func=mod)
jacobian.full(y=c(1,2,3,4),func=mod)


# Boundary value problem - see also the vignette
bvp22 <- function (y, xi) {
  dy2 <- diff(diff(c(ya,y,yb))/dx)/dx
  dy  <- 0.5*(diff(c(ya,y))/dx+diff(c(y,yb))/dx)

  return(xi*dy2+dy+y^2)
}

dx <- 0.001
x  <- seq(0,1,by=dx)
N <- length(x)
ya <- 0
yb <- 0.5

print(system.time(
Y1<- multiroot.1D(f=bvp22, start=runif(N), nspec=1, xi=0.1)
)*1000)
Y2<- multiroot.1D(f=bvp22, start=runif(N), nspec=1, xi=0.05)
Y3<- multiroot.1D(f=bvp22, start=runif(N), nspec=1, xi=0.01)

plot(x,Y3$root, type="l", col="green", lwd=2,
  main="bvp test problem 22" ,ylab="y")
lines(x,Y2$root, col="red", lwd=2)
lines(x,Y1$root, col="blue", lwd=2)
