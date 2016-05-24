library(mbbefd)
library(fitdistrplus)

#____________________________________________________________
#gbeta
n <- 1e3
nboot <- 100
nboot <- 10
set.seed(12345)
x <- rgbeta(n, 2, 2, 5/2)

initpar <- list(shape0=2, shape1=2, shape2=5/2)

if(FALSE)
{
#____________________________________________________________
#test all computation methods
ctr <- list(trace=0, REPORT=1, maxit=1000)
reslist <- NULL
for(meth in c("BFGS", "Nelder", "CG")) #CG with FR update
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta", optim.method=meth, 
                control=ctr, start=initpar))[3]
  reslist <- c(reslist, list(nograd))
}
for(type in 2:3) #CG with PR or BS updates
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta", optim.method="CG", 
                control=c(ctr, type=type), start=initpar))[3] 
  reslist <- c(reslist, list(nograd)) 
}
fullname <- c("BFGS", "NM", paste("CG", c("FR", "PR", "BS")))
names(reslist) <- fullname

dgbeta2 <- function(x, shape0, shape1, shape2, log=FALSE)
  dgbeta(x, exp(shape0), exp(shape1), exp(shape2), log=log)

initpar2 <- lapply(initpar, log)

for(meth in c("BFGS", "Nelder", "CG")) #CG with FR update
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta2", optim.method="BFGS", 
                          control=ctr, start=initpar2))[3]
  nograd$estimate <- exp(nograd$estimate)
  reslist <- c(reslist, list(nograd))
}
for(type in 2:3) #CG with PR or BS updates
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta2", optim.method="CG", 
                        control=c(ctr, type=type), start=initpar2))[3] 
  nograd$estimate <- exp(nograd$estimate)
  reslist <- c(reslist, list(nograd)) 
}
names(reslist)[(length(fullname)+1):length(reslist)] <- paste("exp.", fullname)

getval <- function(x)
  c(x$estimate, loglik=x$loglik, x$counts, x$time)

resNM <- sapply(reslist[grep("NM", names(reslist))], getval)
resCG <- sapply(reslist[grep("CG", names(reslist))], getval)
resBFGS <- sapply(reslist[grep("BFGS", names(reslist))], getval)  
rownames(resNM) <- rownames(resCG) <- rownames(resBFGS) <- c(paste("fitted", c("shape0", "shape1", "shape2")), "fitted loglik", "func. eval. nb.", "grad. eval. nb.", "time (sec)")


#____________________________________________________________
#empirical check of the log-likelihood computation
head(cbind(do.call("dgbeta", c(list(x), initpar, NULL, log=TRUE)), 
log(do.call("dgbeta", c(list(x), initpar, NULL) ) )
))
colSums(cbind(do.call("dgbeta", c(list(x), initpar, NULL, log=TRUE)), 
              log(do.call("dgbeta", c(list(x), initpar, NULL) ) )
))

#____________________________________________________________
#test with starting values equal theoretical values

f1 <- fitDR(x, "oigbeta", "mle")

f1 <- fitdist(x, "gbeta", method="mle", start=initpar) # , control=list(trace=3, REPORT=1)) 
summary(f1)
cdfcomp(f1, do.points=FALSE, ylogscale = TRUE)
lines(0:100/100, pgbeta(0:100/100, 2, 2, 5/2), col="green")

denscomp(f1)
lines(0:100/100, dgbeta(0:100/100, 2, 2, 5/2), col="green")

#____________________________________________________________
#look at the log-likelihood function around the estimated value
par(mfrow=c(1,3))
llsurface(plot.min=c(0.1, 0.1), plot.max=c(5, 4), nlevels=20,
                         plot.arg=c("shape1", "shape2"), fix.arg=as.list(f1$estimate[1]),
                         plot.np=50, obs=x, distr="gbeta", plot.type="contour")
points(f1$estimate["shape1"], f1$estimate["shape2"], pch="+", col="red")
points(2, 5/2, pch="x", col="green")
llsurface(plot.min=c(0.1, 0.1), plot.max=c(6, 4), nlevels=20,
       plot.arg=c("shape0", "shape2"), fix.arg=as.list(f1$estimate[2]),
       plot.np=50, obs=x, distr="gbeta", plot.type="contour")
points(f1$estimate["shape0"], f1$estimate["shape2"], pch="+", col="red")
points(2, 5/2, pch="x", col="green")
llsurface(plot.min=c(0.1, 0.1), plot.max=c(5, 6), nlevels=20,
       plot.arg=c("shape1", "shape0"), fix.arg=as.list(f1$estimate[3]),
       plot.np=50, obs=x, distr="gbeta", plot.type="contour")
points(f1$estimate["shape1"], f1$estimate["shape0"], pch="+", col="red")
points(2, 2, pch="x", col="green")


par(mfrow=c(1,3))
llcurve(plot.min=0.1, plot.max=5, plot.arg="shape0", fix.arg=as.list(f1$estimate[-1]), plot.np=50, 
        obs=x, distr="gbeta", enhance=FALSE)
abline(v=c(2, f1$estimate["shape0"]), col=c("green", "red"))
llcurve(plot.min=0.1, plot.max=4, plot.arg="shape1", fix.arg=as.list(f1$estimate[-2]), plot.np=50, 
                      obs=x, distr="gbeta", enhance=FALSE)
abline(v=c(2, f1$estimate["shape1"]), col=c("green", "red"))
llcurve(plot.min=0.1, plot.max=4, plot.arg="shape2", fix.arg=as.list(f1$estimate[-3]), plot.np=50, 
                      obs=x, distr="gbeta", enhance=FALSE)
abline(v=c(5/2, f1$estimate["shape2"]), col=c("green", "red"))


#bootstrap
b1 <- bootdist(f1, niter=nboot, silent=TRUE)
summary(b1)

plot(b1, enhance=TRUE, trueval=c(2, 2, 5/2))
}

#____________________________________________________________
#init value




s00 <- optimize(function(z) 
  (mbbefd:::Theil.emp(x) - mbbefd:::Theil.theo.shape0(z, obs=x))^2, lower=0.01, upper=20)$minimum
initpar1 <- c(list(shape0=1), as.list(fitdist(x, "beta", method="mme")$estimate))
initpar2 <- c(list(shape0=s00), as.list(fitdist(x^s00, "beta", method="mme")$estimate))




#____________________________________________________________

fitdist(x, "gbeta", method="mle", start=initpar1)
fitdist(x, "gbeta", method="mle", start=initpar2)

f2 <- fitdist(x, "gbeta", method="mle", start=initpar1)
summary(f2)
cdfcomp(f2, do.points=FALSE, ylogscale = TRUE)
lines(0:100/100, pgbeta(0:100/100, 2, 2, 5/2), col="green")
