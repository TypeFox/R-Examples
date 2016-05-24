library(rockchalk)


set.seed(23452345)
N <- 100
x1 <- gl(12, 2, labels = LETTERS[1:12])
x2 <- gl(8, 3, labels = LETTERS[12:24])
x1 <- sample(x = x1, size=N, replace = TRUE)
x2 <- sample(x = x2, size=N, replace = TRUE)
z1 <- rnorm(N)
a1 <- rnorm(N, mean = 1.2, sd = 1.7)
a2 <- rpois(N, lambda = 10 + a1)
a3 <- rgamma(N, 0.5, 4)
b1 <- rnorm(N, mean = 1.3, sd = 1.4)
dat <- data.frame(z1, a1, x2, a2, x1, a3, b1)
summary(dat)


summarize(dat)


summarizeNumerics(dat)
summarizeFactors(dat, maxLevels = 5)

summarize(dat, alphaSort = FALSE)

summarize(dat, digits = 6, alphaSort = FALSE)

summarize(dat, digits = 22, alphaSort = FALSE)

summarize(dat, maxLevels = 2)

datsumm <- summarize(dat)

datsumm$numerics
datsumm[[1]]  ## same: gets numerics

datsumm$factors
datsumm[[2]]


## Use numerics output to make plots. First,
## transpose gives varnames x summary stat matrix
datsummNT <- t(datsumm$numerics)
datsummNT <- as.data.frame(datsummNT)

plot(datsummNT$mean, datsummNT$var, xlab = "The Means",
    ylab = "The Variances")

plot(datsummNT$mean, datsummNT$var, xlab = "The Means",
    ylab = "The Variances", type = "n")
text(datsummNT$mean, datsummNT$var, labels = rownames(datsummNT))

## Here's a little plot wrinkle.  Note variable names are 'out to the
##  edge' of the plot. If names are longer they don't stay inside
##  figure. See?

## Make the variable names longer

rownames(datsummNT)
rownames(datsummNT) <- c("boring var", "var with long name",
    "tedious name var", "stupid varname", "buffoon not baboon")
plot(datsummNT$mean, datsummNT$var, xlab = "The Means",
    ylab = "The Variances", type = "n")
text(datsummNT$mean, datsummNT$var, labels = rownames(datsummNT),
    cex = 0.8)
## That's no good. Names across the edges

## We could brute force the names outside the edges like
##  this
par(xpd = TRUE)
text(datsummNT$mean, datsummNT$var, labels = rownames(datsummNT),
    cex = 0.8)
## but that is not much better
par(xpd = FALSE)

## Here is one fix. Make the unused space inside the plot
##  larger by
## making xlim and ylim bigger.  I use the magRange
##  function from
## rockchalk to easily expand range to 1.2 times its
##  current size.
## otherwise, long variable names do not fit inside plot.
##  magRange
## could be asymmetric if we want, but this use is
##  symmetric.

rownames(datsummNT)
rownames(datsummNT) <- c("boring var", "var with long name",
    "tedious name var", "stupid varname", "buffoon not baboon")
plot(datsummNT$mean, datsummNT$var, xlab = "The Means",
    ylab = "The Variances", type = "n", xlim = magRange(datsummNT$mean,
        1.2), ylim = magRange(datsummNT$var, 1.2))
text(datsummNT$mean, datsummNT$var, labels = rownames(datsummNT),
    cex = 0.8)

## Here's another little plot wrinkle.  If we don't do that to keep
## the names in bounds, we need some fancy footwork.  Note when a
## point is near the edge, I make sure the text prints toward the
## center of the graph.
plot(datsummNT$mean, datsummNT$var, xlab = "The Means",
    ylab = "The Variances")
## calculate label positions. This is not as fancy as it could be.  If
##  there were lots of variables, we'd have to get smarter about
##  positioning labels on above, below, left, or right.
labelPos <- ifelse(datsummNT$mean - mean(datsummNT$mean,
    na.rm = TRUE) > 0, 2, 4)
text(datsummNT$mean, datsummNT$var, labels = rownames(datsummNT),
    cex = 0.8, pos = labelPos)



x <- data.frame(x = rnorm(N), y = gl(50, 2), z = rep(1:4,
    25), ab = gl(2, 50))

summarize(x)
summarize(x, maxLevels = 15)

sumry <- summarize(x)
sumry[[1]]  ##another way to get the numerics output
sumry[[2]]  ##another way to get the factors output

dat <- data.frame(x = rnorm(N), y = gl(50, 2), z = factor(rep(1:4,
    25), labels = c("A", "B", "C", "D")), animal = factor(ifelse(runif(N) <
    0.2, "cow", ifelse(runif(N) < 0.5, "pig", "duck"))))

summarize(dat)

## Run this if you have internet access

## fn <- "http://pj.freefaculty.org/guides/stat/DataSets/USNewsCollege/USNewsCollege.csv"
## dat <- read.table(url(fn), sep = ",")

## colnames(dat) <- c("fice", "name", "state", "private", "avemath",
##                    "aveverb", "avecomb", "aveact", "fstmath",
##                    "trdmath", "fstverb", "trdverb", "fstact",
##                    "trdact", "numapps", "numacc", "numenr",
##                    "pctten", "pctquart", "numfull", "numpart",
##                    "instate", "outstate", "rmbrdcst", "roomcst",
##                    "brdcst", "addfees", "bookcst", "prsnl",
##                    "pctphd", "pctterm", "stdtofac", "pctdonat",
##                    "instcst", "gradrate")

## dat$private <- factor(dat$private, labels = c("public",
##                                    "private"))
## sumry <- summarize(dat, digits = 2)
## sumry

## sumry[[1]]
## sumry[[2]]

## summarize(dat[, c("fice", "name", "private", "fstverb",
##                   "avemath")], digits = 4)

