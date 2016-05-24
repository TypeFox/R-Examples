library(diptest)

stopifnot(dip(c(1,1,2,2)) == 1/4)# the maximal value possible: two point dist

## very first small "unimodal" example --- the  1/(2*n) result:
n <- length(u <- cumsum(0:3))
d <- dip(u, debug=TRUE)# shows the final if() {added by MM} is really needed
stopifnot(d == dip(-u), d == 1/(2*n))# exact "=" for n = 4 !
## Note that I believe this should *not* give 0 (as fmechler@.. did),
## but rather 1/(2n) because that's  (1/n) / 2  and
## (1/n) is the correct distance between LCM and GCM

## Small example -- but MM sees difference (32-bit / 64-bit):
x <- c(0,2:3,5:6)
d1 <- dip(x,   full=TRUE, debug=2)
d2 <- dip(6-x, full=TRUE, debug=2)
str(d1)
str(d2)

if(!dev.interactive(orNone=TRUE)) pdf("ex1.pdf")
par(mfrow = 2:1, mar = .1+c(3,4,2,1), mgp=c(1.5,.6,0), oma = c(0,0,2.1,0))
#
plot(d1)
abline(v=-1:7, h = seq(0,1,by=0.2), lty="83", col = "gray")
#
plot(d2)
abline(v=-1:7, h = seq(0,1,by=0.2), lty="83", col = "gray")
#
## "title" only now
mtext("dip() problem with 'mirror x'", side=3, line = 0.8,
      outer=TRUE, cex = 1.5, font = 2)


##  Yong Lu <lyongu+@cs.cmu.edu> example -- a bit smaller
x2 <- c(1, rep(2, 9))
stopifnot(dip(x2) == dip(3 - x2))
str(dip(x2, full=TRUE))
cat('Time elapsed: ', (.pt <- proc.time()),'\n') # "stats"

## Real data examples :

data(statfaculty)

str(dip(statfaculty, full = "all", debug = 3), vec.len = 8)

data(faithful)
fE <- faithful$eruptions
str(dip(fE, full = "all", debug = 3),
    vec.len= 8)

data(precip)
str(dip(precip, full = TRUE, debug = TRUE))

cat('Time elapsed: ', proc.time() - .pt,'\n') # "stats"

if(!interactive()) warnings()
