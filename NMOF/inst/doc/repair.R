### R code from vignette source 'repair.Rnw'

###################################################
### code chunk number 1: repair.Rnw:24-25
###################################################
options(continue = " ", digits = 5)


###################################################
### code chunk number 2: repair.Rnw:42-44
###################################################
set.seed(112233)
options(digits = 3)


###################################################
### code chunk number 3: repair.Rnw:58-62
###################################################
up <- rep(1, 4L)
lo <- rep(0, 4L)
x <- rnorm(4L)
x


###################################################
### code chunk number 4: repair.Rnw:65-69
###################################################
repair1a <- function(x,up,lo) 
    pmin(up,pmax(lo,x))
x
repair1a(x, up, lo)


###################################################
### code chunk number 5: repair.Rnw:74-88
###################################################
repair1b <- function(x, up, lo) {
    ii <- x > up
    x[ii] <- up[ii]
    ii <- x < lo
    x[ii] <- lo[ii]
    x
}
repair1c <- function(x, up, lo) {
    xadjU <- x - up
    xadjU <- xadjU + abs(xadjU)
    xadjL <- lo - x
    xadjL <- xadjL + abs(xadjL)
    x - (xadjU - xadjL)/2
}


###################################################
### code chunk number 6: repair.Rnw:98-107
###################################################
repair1a(x, up, lo)
repair1b(x, up, lo)
repair1c(x, up, lo)

trials <- 5000L
strials <- seq_len(trials)
system.time(for(i in strials) y1 <- repair1a(x, up, lo))
system.time(for(i in strials) y2 <- repair1b(x, up, lo))
system.time(for(i in strials) y3 <- repair1c(x, up, lo))


###################################################
### code chunk number 7: repair.Rnw:111-114
###################################################
X <- array(rnorm(25L), dim = c(5L, 5L))
X
repair1c(X, up = 0.5, lo = -0.5)


###################################################
### code chunk number 8: repair.Rnw:120-124
###################################################
pmax2 <- function(x1, x2) 
    ((x1 + x2) + abs(x1 - x2)) / 2
pmin2 <- function(x1, x2) 
    ((x1 + x2) - abs(x1 - x2)) / 2


###################################################
### code chunk number 9: repair.Rnw:127-139
###################################################
x1 <- rnorm(100L)
x2 <- rnorm(100L)

t1 <- system.time(for (i in strials) z1 <- pmax(x1,x2) )
t2 <- system.time(for (i in strials) z2 <- pmax2(x1,x2))
t1[[3L]]/t2[[3L]] ## speedup
all.equal(z1, z2)

t1 <- system.time(for (i in strials) z1 <- pmin(x1,x2) )
t2 <- system.time(for (i in strials) z2 <- pmin2(x1,x2))
t1[[3L]]/t2[[3L]] ## speedup
all.equal(z1, z2)


###################################################
### code chunk number 10: repair.Rnw:149-174
###################################################
repair2 <- function(x, up, lo) {
    done <- TRUE
    e <- sum(x - up  + abs(x - up) + lo - x  + abs(lo - x))
    if (e > 1e-12) 
        done <- FALSE
    r <- up - lo
    while (!done) {
        adjU <- x - up
        adjU <- adjU + abs(adjU)
        adjU <- adjU + r - abs(adjU - r)

        adjL <- lo - x
        adjL <- adjL + abs(adjL)
        adjL <- adjL + r - abs(adjL - r)

        x <- x - (adjU - adjL)/2
        e <- sum(x - up  + abs(x - up) + lo - x  + abs(lo - x))
        if (e < 1e-12) 
            done <- TRUE
    }
    x
}
x
repair2(x, up, lo)
system.time(for (i in strials) y4 <- repair2(x,up,lo))


###################################################
### code chunk number 11: repair.Rnw:194-198
###################################################
T <- 20L
x <- logical(T)
x[runif(T) < 0.4] <- TRUE
x


###################################################
### code chunk number 12: repair.Rnw:202-204
###################################################
kmax <- 5L
kmin <- 3L


###################################################
### code chunk number 13: repair.Rnw:208-223
###################################################
resample <- function(x, ...) x[sample.int(length(x), ...)]
repairK <- function(x, kmax, kmin) {
    sx <- sum(x)
    if (sx > kmax) {
        i <- resample(which(x), sx - kmax)
        x[i] <- FALSE
    } else if (sx < kmin) {
        i <- resample(which(!x), kmin - sx)
        x[i] <- TRUE
    }
    x
}
printK <- function(x)
    cat(paste(ifelse(x, "o", "."), collapse = ""),
        "-- cardinality", sum(x), "\n")


###################################################
### code chunk number 14: repair.Rnw:226-231
###################################################
for (i in 1:10) {
    if (i==1L) printK(x)
    x1 <- repairK(x, kmax, kmin)
    printK(x1)
}


###################################################
### code chunk number 15: repair.Rnw:234-240
###################################################
x <- logical(T); x[10L] <- TRUE
for (i in 1:10) {
    if (i==1L) printK(x)
    x1 <- repairK(x, kmax, kmin)
    printK(x1)
}


