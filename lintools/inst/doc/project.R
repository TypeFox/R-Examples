## ----echo=FALSE,fig.width=5,fig.height=5---------------------------------
x <- c(-0.5,2)
y1 <- x
y2 <- 1-x
plot( x,y1,'l',xlim=c(-0.5,1.5),ylim=c(-0.5,1.5),xlab='x',ylab='y',lwd=1.8)
lines(x,y2,lwd=2)
polygon(
    x = c(-0.5, 1.5,-.5)
  , y = c(-0.5, 1.5,1.5)
  , density = 5
  , border=NA
  , angle=-45
)

polygon(
  x = c(-0.5,1.5,1.5)
  ,y =c(1.5,-0.5,1.5) 
  ,density=5
  ,border=NA
)
abline(h=0,v=0)

points(0.8,-0.2,pch=16,cex=1.8)
arrows(x0=0.8,y=-0.2,x1=1,y1=0,length=0.1,lwd=2)
points(0.5,0.5,pch=16,col='grey',cex=1.8)
arrows(x0=1,y0=0,x1=0.5,y1=0.5,length=0.1,lwd=2)

text(-0.15,-0.3,labels="y = x")
text(1.05,-0.3,labels="y = 1 - x")

## ------------------------------------------------------------------------
library(lintools)
x <- c(0.8,-0.2)
A <- matrix(c(1,-1,-1,-1), byrow=TRUE, nrow=2)
b <- c(0,-1)

## ------------------------------------------------------------------------
project(x=x,A=A,b=b,neq=0)

## ------------------------------------------------------------------------
A <- data.frame(
  row = c(1,1,2,2)
  ,col = c(1,2,1,1)
  ,coef = c(1,-1,-1,-1)
)
b <- c(0,-1)
x <- c(0.8,-0.2)

## ------------------------------------------------------------------------
sparse_project(x, A=A, b=b, neq=0)

## ------------------------------------------------------------------------
sc <- sparse_constraints(A,b,neq=0)

## ------------------------------------------------------------------------
sc$project(x=c(0.8,-0.2))
# the same problem, but with differing weights
sc$project(x=c(0.8,-0.2),w=c(1,10))

