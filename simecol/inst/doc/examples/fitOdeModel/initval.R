## =======================================================
## Example how to fit initial values
## =======================================================

library(simecol)

## ======== load example model =========
data(chemostat)

## derive scenario
cs1 <- chemostat

## generate some noisy data
set.seed(12571) # make example reproducible
parms(cs1)[c("vm", "km")] <- c(1, 5)
times(cs1) <- c(from=0, to=40, by=2)
obs <- out(sim(cs1))

## add some noise to the simulated data
obs$S <- obs$S + rnorm(obs$S, sd= 0.2 * sd(obs$S))
obs$X <- obs$X + rnorm(obs$X, sd= 0.1 * sd(obs$X))

## split obs to yobs and time
obstime <- obs$time
yobs <- obs[c("X", "S")]

## add initial values to the parameter vector
parms(cs1) <- c(parms(cs1), X=10, S=10)

## define an intifunc that copies these parameters back to init
initfunc(cs1) <- function(obj) {
  init(obj) <- parms(obj)[c("X", "S")] # Note!  Order is important!
  obj
}

## set *external* time step to same as in observations,
## and use efficient algorithm with automatic *internal* time steps
times(cs1) <- obstime
solver(cs1) <- "lsoda"

whichpar  <- c("vm", "km", "D", "X")
parms(cs1)[whichpar] <- c(vm=1, km=1, D=1, X=5)

lower <- c(vm=0.1, km=0.1, D=0, X=0.1)
upper <- c(vm=10, km=40, D=5, X=200)


res.p <- fitOdeModel(cs1, whichpar = whichpar, obstime, yobs,
                     debuglevel=0, fn = ssqOdeModel,
                     method = "PORT", lower = lower, upper = upper,
                     #control=list(trace=TRUE),
                     atol=1e-4, rtol=1e-4)

res.m <- fitOdeModel(cs1, whichpar = whichpar, obstime, yobs,
                     debuglevel=0, fn = ssqOdeModel,
                     method = "Nelder-Mead", lower = lower, upper = upper,
                     #control=list(trace=TRUE),
                     atol=1e-4, rtol=1e-4)


res.b <- fitOdeModel(cs1, whichpar = whichpar, obstime, yobs,
                     debuglevel=0, fn = ssqOdeModel,
                     method = "bobyqa", #lower = lower, upper = upper,
                     #control=list(iprint=2),
                     atol=1e-4, rtol=1e-4)

res.n <- fitOdeModel(cs1, whichpar = whichpar, obstime, yobs,
                     debuglevel=0, fn = ssqOdeModel,
                     method = "newuoa", #lower = lower, upper = upper,
                     #control=list(iprint=2),
                     atol=1e-4, rtol=1e-4)



## set small external time step for good graphics
times(cs1) <- seq(0, 40, length.out=200)

## assign fitted parameters to scenarios
cs.b <- cs.n <- cs.p <- cs.m <- cs1

parms(cs.b)[whichpar] <- coef(res.b)
parms(cs.p)[whichpar] <- coef(res.p)
parms(cs.m)[whichpar] <- coef(res.m)
parms(cs.n)[whichpar] <- coef(res.n)

cs.p <- sim(cs.p)
cs.m <- sim(cs.m)
cs.n <- sim(cs.n)
cs.b <- sim(cs.b)


## compare results
plot(cs.p, cs.m, cs.n, cs.b, obs=obs)
legend("bottomright", legend=c("Port", "Nelder-Mead", "newuoa", "bobyqa"), lty=1:4, col=1:4)
