# Following plots and demos illustrate various functionality from the package.
# Some of the parameters have been altered in the interest of speed.


# Change 'n' to higher values to see higher resolution plots.
n <- 200

# Change 'val' to see the functions evaluated over a larger complex area
val <- 4

# Change 'no.of.schemes' to see more schemes used
no.of.schemes <- 3

# Change 'fact' to adjust the factor used for the upper-half-plane functions
fact <- 3


# Remember the settings:
opar <- par(ask = dev.interactive(orNone = TRUE))



# Define a complex grid:
x <- seq(from = -val, to=val, len=n)
y <- x
z <- outer(x,1i*x,"+")

# Now a grid for the upper half-plane functions:
xupper <- x/fact
yupper <- (y+val+val/n)/fact
zupper <- outer(xupper,1i*yupper,"+")


# A little wrapper for view():
f <- function(...){
  view(... , axes=FALSE,xlab="",ylab="")
  jj <- c(-val,-val/2,0,val/2,val)
  axis(1,pos=-val,at=jj)
  axis(2,pos=-val,at=jj)
}

# A wrapper for view() for the upper half-plane functions:
fupper <- function(...){
  view(... , axes=FALSE,xlab="",ylab="")
  jj <- c(-val,val)/fact
  jj2 <- c(0,val*2)/fact
  axis(1,pos=0,at=pretty(jj))
  axis(2,pos=-val/fact,at=pretty(jj2))
}

# Tiny function for the title:
fish <- function(string,i){
  paste(string, ". scheme=", as.character(i),"; try increasing n",sep="")
}

# Wrapper to run view() a few times with differing args:
jj <- function(fz,string){
  for(i in sample(19,no.of.schemes)){
    f(x,y,fz,scheme=i,real=FALSE,main=fish(string,i))
  }
  f(x,y,fz,scheme=-1,real=TRUE,imag=FALSE,nlevels=33,drawlabels=FALSE,main=fish(string,-1))
  f(x,y,fz,scheme=-1,real=TRUE,nlevels=33,drawlabels=FALSE,main=fish(string,-1))
}

# Corresponding wrapper for the upper half plane functions:
kk <- function(fz,string){
  for(i in sample(19,no.of.schemes)){
    fupper(xupper,yupper,fz,scheme=i,real=FALSE,main=fish(string,i))
  }
  fupper(xupper,yupper,fz,scheme=-1,real=TRUE,imag=FALSE,
         nlevels=33,drawlabels=FALSE,main=fish(string,-1))
  fupper(xupper,yupper,fz,scheme=-1,real=TRUE,imag=TRUE,
         nlevels=33,drawlabels=FALSE,main=fish(string,-1))
}

# Now run everything; jj() and kk() take some time to complete:
jj(fz=limit(sn(z,m=1/2+0.6i)),"sn(z)")
jj(fz=limit(P(z,c(1+2.1i,1.3-3.2i))),"P(z)")
jj(fz=limit(zeta(z,c(1+1i,2-3i))),"zeta(z)")
jj(fz=limit(sigma(z,c(10+11i,20-31i))),"sigma(z)")
kk(fz=limit(J(zupper,maxiter=100)),"J(z)")
kk(fz=limit(lambda(zupper,maxiter=100)),"lambda(z)")



# reset old settings:
par(opar)
