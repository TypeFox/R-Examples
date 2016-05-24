## M: computer model test functon used in Goh et al, 2013 (Technometrics)
## an elaboration of one from Bastos and O'Hagan, 2009 (Technometrics) 
M <- function(x,u) 
  {
    x <- as.matrix(x)
    u <- as.matrix(u)
    out <- (1-exp(-1/(2*x[,2]))) 
    out <- out * (1000*u[,1]*x[,1]^3+1900*x[,1]^2+2092*x[,1]+60) 
    out <- out / (100*u[,2]*x[,1]^3+500*x[,1]^2+4*x[,1]+20)  
    return(out)
  }
  
## bias: discrepancy function from Goh et al, 2013 
bias <- function(x) 
  {
    x<-as.matrix(x)   
    out<- 2*(10*x[,1]^2+4*x[,2]^2) / (50*x[,1]*x[,2]+10)
    return(out)
  }

## beta.prior: marginal beta prior for u, 
## defaults to a mode at 1/2 in hypercube
beta.prior <- function(u, a=2, b=2, log=TRUE)
{
  if(length(a) == 1) a <- rep(a, length(u))
  else if(length(a) != length(u)) stop("length(a) must be 1 or length(u)")
  if(length(b) == 1) b <- rep(b, length(u))
  else if(length(b) != length(u)) stop("length(b) must be 1 or length(u)")
  if(log) return(sum(dbeta(u, a, b, log=TRUE)))
  else return(prod(dbeta(u, a, b, log=FALSE)))
}


## tgp for LHS sampling
library(tgp)
rect <- matrix(rep(0:1, 4), ncol=2, byrow=2)

## training and testing inputs
ny <- 50; nny <- 1000  
X <- lhs(ny, rect[1:2,])    ## computer model train
XX <- lhs(nny, rect[1:2,],) ## test

## true (but unknown) setting of the calibration parameter
## for the computer model
u <- c(0.2, 0.1)
Zu <- M(X, matrix(u, nrow=1)) 
ZZu <- M(XX, matrix(u, nrow=1)) 

## field data response, biased and replicated
sd <- 0.5
## Y <- computer output + bias + noise
reps <- 2 ## example from paper uses reps <- 10
Y <- rep(Zu,reps) + rep(bias(X),reps) + rnorm(reps*length(Zu), sd=sd) 
YYtrue <- ZZu + bias(XX) 
## variations: remove the bias or change 2 to 1 to remove replicates

## computer model design
nz <- 10000
XU <- lhs(nz, rect)
nth <- 8 ## number of threads to use in emulation

## augment with physical model design points 
## with various u settings
XU2 <- matrix(NA, nrow=10*ny, ncol=4)
for(i in 1:10) {
  I <- ((i-1)*ny+1):(ny*i)
  XU2[I,1:2] <- X
}
XU2[,3:4] <- lhs(10*ny, rect[3:4,])
XU <- rbind(XU, XU2)

## evaluate the computer model
Z <- M(XU[,1:2], XU[,3:4])

## flag indicating if estimating bias/discrepency or not
bias.est <- TRUE
## two passes of ALC with MLE calculations for aGP.seq
methods <- rep("alc", 2)

## set up priors
da <- d <- darg(NULL, XU)
g <- garg(list(mle=TRUE), Y) 

## initial space-filling search for starting value of
## the calibration parameter, u
initsize <- 10*2 
imesh <- 0.1
irect <- rect[1:2,]
irect[,1] <- irect[,1] + imesh/2
irect[,2] <- irect[,2] - imesh/2
uinit.cand <- lhs(10*initsize, irect) 
uinit <- dopt.gp(initsize, Xcand=lhs(10*initsize, irect))$XX
llinit <- rep(NA, nrow(uinit))
for(i in 1:nrow(uinit)) {
  llinit[i] <- fcalib(uinit[i,], XU, Z, X, Y, da, d, g, beta.prior, 
                  methods, M, bias.est, nth)
}

## now for a NOMAD search via snomadr in the crs package
library(crs)
opts <- list("MAX_BB_EVAL"=1000, "INITIAL_MESH_SIZE"=imesh, "MIN_POLL_SIZE"="r0.001")
its <- 0
o <- order(llinit)
i <- 1
out <- NULL
while(its < 10) {
  cat("NOMAD start=", i, ", uinit=(", paste(uinit[o[i],], collapse=", "), ")\n", sep="")
  outi <- snomadr(fcalib, 2, c(0,0), 0, x0=uinit[o[i],],
            lb=c(0,0), ub=c(1,1), opts=opts, XU=XU, 
            Z=Z, X=X, Y=Y, da=da, d=d, g=g, methods=methods, M=M, 
            bias=bias.est, omp.threads=nth, uprior=beta.prior, 
            save.global=TRUE, verb=1)
  its <- its + outi$iterations
  ## keep the best
  if(is.null(out) || outi$objective < out$objective) out <- outi
  i <- i + 1;
}

## visualize the posterior surface
library(akima)
Xp <- rbind(uinit, as.matrix(fcalib.save[,1:2]))
Zp <- c(-llinit, fcalib.save[,3])
wi <- which(!is.finite(Zp))
if(length(wi) > 0) { Xp <- Xp[-wi,]; Zp <- Zp[-wi]}
surf <- interp(Xp[,1], Xp[,2], Zp, duplicate="mean")
image(surf, xlab="u1", ylab="u2", main="posterior surface",
  col=heat.colors(128), xlim=c(0,1), ylim=c(0,1))
points(uinit)
points(fcalib.save[,1:2], col=3, pch=18)

## the maximizing u-value does not agree with the true value
u.hat <- outi$solution
points(u.hat[1], u.hat[2], col=4, pch=18)
print(rbind(u.hat, u))
## put u on plot as cross-hairs
abline(v=u[2], lty=2)
abline(h=u[1], lty=2)

## but which one has a higher value of the log posterior?
Xu <- cbind(X, matrix(rep(u, ny), ncol=2, byrow=TRUE))
Mhat.u <- aGP.seq(XU, Z, Xu, da, methods, ncalib=2, omp.threads=nth, verb=0)
cmle.u <- discrep.est(X, Y, Mhat.u$mean, d, g, bias.est, FALSE)
cmle.u$ll <- cmle.u$ll + beta.prior(u)
data.frame(u.hat=-outi$objective, u=cmle.u$ll)
## u.hat has the higher "llik" value

## now lets do a predictive comparison

## first with the true u value
XXu <- cbind(XX, matrix(rep(u, nny), ncol=2, byrow=TRUE))
Mhat.oos.u <- aGP.seq(XU, Z, XXu, da, methods, ncalib=2, 
  omp.threads=nth, verb=0)
YYm.pred.u <- predGP(cmle.u$gp, XX)
YY.pred.u <- YYm.pred.u$mean + Mhat.oos.u$mean
rmse.u <- sqrt(mean((YY.pred.u - YYtrue)^2))
deleteGP(cmle.u$gp)

## now with the estimated u value
Xu <- cbind(X, matrix(rep(u.hat, ny), ncol=2, byrow=TRUE))
Mhat <- aGP.seq(XU, Z, Xu, da, methods, ncalib=2, omp.threads=nth, verb=0)
cmle <- discrep.est(X, Y, Mhat$mean, d, g, bias.est, FALSE)
cmle$ll <- cmle$ll + beta.prior(u.hat)
## should give the same number
print(c(cmle$ll, -outi$objective))
## back to predicting out of sample
XXu <- cbind(XX, matrix(rep(u.hat, nny), ncol=2, byrow=TRUE))
Mhat.oos <- aGP.seq(XU, Z, XXu, da, methods, ncalib=2, omp.threads=nth, verb=0)
YYm.pred <- predGP(cmle$gp, XX)
YY.pred <- YYm.pred$mean + Mhat.oos$mean
rmse <- sqrt(mean((YY.pred - YYtrue)^2))

## compare rmses, u.hat is usually better
data.frame(u.hat=rmse, u=rmse.u)