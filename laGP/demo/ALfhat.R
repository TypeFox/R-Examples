## augmenting the example in ?optim.auglag, here is a
## a variation with an unknown objective, where the objective
## is taken from other laGP experiece, see e.g., ?laGP

## since the demo is stochastic, we suggest re-running it 
## several times to get a sense of typical behavior, especially
## with regards to comparing the statistical surrogate model
## performance (optim.auglag) with a traditional AL approach
## via optim


## set bounding rectangle for adaptive sampling
B <- matrix(c(rep(0,2),rep(1,2)),ncol=2)

## objective function
f2d <- function(x, y=NULL)
  {
    if(is.null(y)) {
      if(!is.matrix(x)) x <- matrix(x, ncol=2)
        y <- x[,2]; x <- x[,1]
    }
    g <- function(z)
        return(exp(-(z-1)^2) + exp(-0.8*(z+1)^2) - 0.05*sin(8*(z+0.1)))
    z <- -g(x)*g(y)
    return(z)
  }

## new full objective-constraints function modified to use f2d
aimprob2 <- function(X, known.only = FALSE)
{
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  if(known.only) stop("no outputs are treated as known")
  f <- f2d(4*(X-0.5))
  c1 <- 1.5 - X[,1] - 2*X[,2] - 0.5*sin(2*pi*(X[,1]^2 - 2*X[,2]))
  c2 <- rowSums(X^2)-1.5
  return(list(obj=f, c=cbind(c1,c2)))
}


## for visualization
x <- seq(0,1, length=200)
X <- expand.grid(x, x)
out <- aimprob2(as.matrix(X))
fv <- out$obj
fv[out$c[,1] > 0 | out$c[,2] > 0] <- NA
fi <- out$obj
fi[!(out$c[,1] > 0 | out$c[,2] > 0)] <- NA
plot(0, 0, type="n", xlim=B[1,], ylim=B[2,], xlab="x1", ylab="x2",
  main="optim.auglag: blue/diamonds; optim: black/triangles")
contour(x, x, matrix(out$c[,1], ncol=length(x)), nlevels=1, levels=0, 
  drawlabels=FALSE, add=TRUE, lwd=2)
contour(x, x, matrix(out$c[,2], ncol=length(x)), nlevels=1, levels=0, 
  drawlabels=FALSE, add=TRUE, lwd=2)
contour(x, x, matrix(fv, ncol=length(x)), nlevels=10, add=TRUE, col="forestgreen")
contour(x, x, matrix(fi, ncol=length(x)), nlevels=13, add=TRUE, col=2, lty=2)

## optimize under constraints
out2 <- optim.auglag(aimprob2, B, fhat=TRUE, ab=c(3/2,8), start=20, end=100)

## put the answer found on the plot
v <- apply(out2$C, 1, function(x) { all(x <= 0) })
X <- out2$X[v,]
obj <- out2$obj[v]
xbest <- X[which.min(obj),]

## Optionally, wrap up by drilling down with a standard AL approach

## augmented Lagrangian objective function
aimprob2.AL <- function(x, B, lambda, rho)
{
  if(any(x < B[,1]) | any(x > B[,2])) return(Inf)
  fc <- aimprob2(x)
  al <- fc$obj + lambda %*% drop(fc$c) + rep(1/(2*rho), 2) %*% pmax(0, drop(fc$c))^2
  return(al) 
}

## loop over AL updates until a valid solution is found
lambda <- out2$lambda[nrow(out2$lambda),]; rho <- out2$rho[length(out2$rho)]
while(1) {
    o <- optim(xbest, aimprob2.AL, control=list(maxit=15),
      B=B, lambda=lambda, rho=rho)
    fc <- aimprob2(o$par)
    points(o$par[1], o$par[2], pch=18, col="blue")
    segments(xbest[1], xbest[2], o$par[1], o$par[2])
    if(all(fc$c <= 0)) { break
    } else {
      lambda <- pmax(0, lambda + (1/rho)*fc$c)
      rho <- rho/2
      xbest <- o$par
    }
  }


## finally, lets see how this compares to a standard AL
## approach initialized randomly, and randomly restarted
## until a total budget is met

lambda <- rep(0,2); rho <- 1/2
xbest <- runif(2)

for(i in 1:6) {
  o <- optim(xbest, aimprob2.AL, control=list(maxit=18),
    B=B, lambda=lambda, rho=rho)
  fc <- aimprob2(o$par)
  points(o$par[1], o$par[2], pch=17)
  segments(xbest[1], xbest[2], o$par[1], o$par[2])

  ## possibly restart
  if(all(fc$c <= 0) && all(xbest == o$par)) {
     xbest <- runif(2); lambda <- rep(0,2); rho <- 1/2
  } else {
    lambda <- pmax(0, lambda + (1/rho)*fc$c)
    if(any(fc$c > 0)) rho <- rho/2
    xbest <- o$par
  }
}