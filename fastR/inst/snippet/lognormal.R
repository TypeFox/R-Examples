#note: these functions do not check for negative values 
#       where they shouldn't be
dlognormal <- function(x,mu=0,sigma=1) {
    dnorm(log(x),mean=mu, sd=sigma) * 1/x
}
rlognormal <- function(n, mu=0, sigma=1) {
    normals <- rnorm(n, mean=mu, sd=sigma )
    return(exp(normals))
}
plognormal <- function(x, mu=0, sigma=1) {
    pnorm( log(x), mean=mu, sd=sigma ) 
}
qlognormal <- function(p, mu=0, sigma=1) {
    exp( qnorm(p, mean=mu, sd=sigma ) )
}
# some checks
randomData <- rlognormal(100,mu=0,sigma=1/2)
quant <- quantile(randomData)
x <- qlognormal(c(0.25,0.5,0.75),mu=0,sigma=1/2); x
plognormal(x,mu=0,sigma=1/2)
plognormal(quant,mu=0,sigma=1/2)

plot1 <- histogram(~randomData)

x <- seq(0,10,by=0.25)
nx <- length(x)
mu <- c(-1, 0, 1)
nmu <- length(mu)
sigma <- c(1/8, 1/4, 1/2, 1, 2, 4)
nsigma <- length(sigma)

x <- rep(x,each=nmu*nsigma)
mu <- rep(rep(mu,nsigma),times=nx)
sigma <- rep(rep(sigma,each=nmu),times=nx)
density <- dlognormal(x,mu,sigma)

plot2 <- xyplot(density~x|paste("sigma", '=',sigma), 
                groups = paste("mu =",mu), 
                type="l", 
                key=simpleKey(paste("mu =",sort(unique(mu))),
                        points=FALSE, lines=TRUE, columns=3),
                scales=list(y=list(relation="free",alternating=FALSE)),
                main = "pdfs of lognormal distributions")
