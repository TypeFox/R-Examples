set.seed(257)
opar <- par(ask = dev.interactive(orNone = TRUE))


cat("******************************************************\n",
    "We give here an example on the mouflon. It just \n",
    "illustrates the use of the functions to handle trajectories\n",
    "with adehabitat",
    "\n******************************************************\n")


##************************************************************
##
## First approach on the data

data(mouflon)
mouflon
id(mouflon)
burst(mouflon)
head(mouflon[[1]])

## plot the data
plot(mouflon)
trajdyn(mouflon)


## Constant time lag between relocations ?
is.regular(mouflon)
mouflon

## yes, one relocation every 20 minutes

## There are some missing values in the data
summaryNAltraj(mouflon)

## But only a weak proportion! Only one for the second week-end!
## Note, however, that for the first week-end, the missing values
## are not randomly distributed:
runsNAltraj(mouflon[1])
plotNAltraj(mouflon[1])

## There is a runs of missing values on sunday... to be kept
## in mind for the rest of this analysis (though the weak
## proportion of missing values is here unlikely to affect the
## results of the analysis)


## For the sake of simplicity, we will work only on the first
## week-end:
mouflon <- mouflon[1]



##************************************************************
##
## Comparison with a model:


## A correlated random walk is made of successive "steps"
## independent from the point of view of the step length
## and the relative angles between successive relocations

## Does this model describes the trajectories adequately?
## we test the independence of relative angles between
## successive steps:

testang.ltraj(mouflon, "relative")

## Not significant

## Are the steps length independent?
wawotest(mouflon)
## There is a positive autocorrelation.
## Another way to perform this test:
indmove(mouflon)


## Time series approach:
## The autocorrelation function on
## the distance:
acf(na.omit(mouflon[[1]]$dist))
## A lack of independence with lag = 1


## Look at the periodogram
spectrum(na.omit(mouflon[[1]]$dist), taper=0, log="no")
## No clear pattern emerges: no periodicity

## A time plot of the step length
plotltr(mouflon, "dist")

## This is not white noise! what already indicated the tests...
## There are some periods when the animal is moving more slowly
## than others.



##************************************************************
##
## Segmentation: STILL UNDER RESEARCH!!!!!!!

## We try to partition the trajectory into several types of behaviours
## we will work on the dist, and suppose a chi distribution for
## the distribution of distances, with different scaling factors

## The models will be chi distribution with different scaling
## factors:
## The function foo allows to estimate the scaling factor for a
## chi distribution from a data frame containing the dx, dy and dt
## component of a trajectory (see hbrown).

foo <- function(x)
{
    u1 <- x$dx/sqrt(x$dt)
    u2 <- x$dy/sqrt(x$dt)
    oo <- cbind(u1,u2)
    oo <- oo[!apply(oo,1,function(y) any(is.na(y))),]
    vc <- crossprod(scale(oo, scale = FALSE))/nrow(oo)
    h <- sqrt(mean(diag(vc)))
    return(h)
}

## Compute the scaling factor for sliding window
sliwinltr(mouflon, foo, step=5, type="locs")

## OK, we have roughly three types of movements:
scfac <- c(0.5, 1.5, 2.5)
(limod <- as.list(paste("dchi(dist/(sqrt(dt)*",
                        scfac, "))")))

## Then, build the probability matrix
mod1 <- modpartltraj(mouflon, limod)
mod1

## Computes the optimal number of segments
bestpartmod(mod1, Km=70)


## The partition with 21 segments
par <- partmod.ltraj(mouflon, 21, mod1)
par
plot(par)

## In the case where the periodogram on the distance indicates a
## period in the movements of the animal, it could be of interest
## to study the periodicity of the types of movements (are
## certain types of movements preferably used at certain periods of
## the day?). Indeed, there is no hypothesis of stationarity
## With this partitioning algorithm

## A possible analysis could then be to identify correlates
## between the type of movements and the habitat



##************************************************************
##
## Rediscretization

## Another common way to analyse trajectory is to rediscretize them
## with a constant step length

plot(mouflon)
red <- redisltraj(mouflon, 100)
plot(red)

## ...Note that the trajectory is no longer regular
red

## We do not present this type of analysis here, as a
## more intensively sampled trajectory is needed for this
## type of analysis


cat("*******************************************************\n",
    "The deeply commented source for this demo can be found in the file:\n",
    file.path(system.file(package = "adehabitat"), "demo", "analysisltraj.r\n"),
    "Examples of management of trajectories are given in demo(managltraj)\n",
    "******************************************************\n")
