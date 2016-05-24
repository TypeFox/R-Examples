# Avoid random failures to converge
set.seed(1)

# See Becker-Clogg (1989) test

library(logmult)
data(color)

rcm <- rc(color[,,1], 2, weighting="marginal", start=NA)
rcu <- rc(color[,,1], 2, weighting="uniform", start=NA)
rcn <- rc(color[,,1], 2, weighting="none", start=NA)

stopifnot(all.equal(fitted(rcm), fitted(rcu)))
stopifnot(all.equal(fitted(rcm), fitted(rcn)))

phim <- maor(fitted(rcm), TRUE, weighting="marginal", norm=2)
phiu <- maor(fitted(rcu), TRUE, weighting="uniform", norm=2)
phin <- maor(fitted(rcn), TRUE, weighting="none", norm=2)

cphim <- maor(fitted(rcm), TRUE, TRUE, weighting="marginal", norm=2)
cphiu <- maor(fitted(rcm), TRUE, TRUE, weighting="uniform", norm=2)
cphin <- maor(fitted(rcm), TRUE, TRUE, weighting="none", norm=2)

maorm <- maor(fitted(rcm), weighting="marginal", norm=2)
maoru <- maor(fitted(rcu), weighting="uniform", norm=2)
maorn <- maor(fitted(rcn), weighting="none", norm=2)

cmaorm <- maor(fitted(rcm), cell=TRUE, weighting="marginal", norm=2)
cmaoru <- maor(fitted(rcu), cell=TRUE, weighting="uniform", norm=2)
cmaorn <- maor(fitted(rcn), cell=TRUE, weighting="none", norm=2)

stopifnot(all.equal(phim, sqrt(sum((rcm$assoc$phi)^2))))
stopifnot(all.equal(phiu, sqrt(sum(abs(rcu$assoc$phi)^2))))
stopifnot(all.equal(phin, sqrt(sum(abs(rcn$assoc$phi)^2))))

stopifnot(all.equal(phim, sqrt(sum(cphim))))
stopifnot(all.equal(phiu, sqrt(sum(cphiu))))
stopifnot(all.equal(phin, sqrt(sum(cphin))))

stopifnot(all.equal(maorm, exp(sqrt(sum(cmaorm)))))
stopifnot(all.equal(maoru, exp(sqrt(sum(cmaoru)))))
stopifnot(all.equal(maorn, exp(sqrt(sum(cmaorn)))))


# Test on perfectly symmetric association
data(ocg1973)

rcm <- rc(ocg1973, 2, symmetric=TRUE, weighting="marginal", start=NA)
rcu <- rc(ocg1973, 2, symmetric=TRUE, weighting="uniform", start=NA)
rcn <- rc(ocg1973, 2, symmetric=TRUE, weighting="none", start=NA)

stopifnot(all.equal(fitted(rcm), fitted(rcu)))
stopifnot(all.equal(fitted(rcm), fitted(rcn)))

w <- (rcm$assoc$row.weights + rcm$assoc$col.weights)[,1]/2
phim <- maor(fitted(rcm), TRUE, weighting="marginal", norm=2, row.weights=w, col.weights=w)
phiu <- maor(fitted(rcu), TRUE, weighting="uniform", norm=2)
phin <- maor(fitted(rcn), TRUE, weighting="none", norm=2)

sphim <- maor(fitted(rcm), TRUE, component="symmetric", weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
sphiu <- maor(fitted(rcu), TRUE, component="symmetric", weighting="uniform", norm=2)
sphin <- maor(fitted(rcn), TRUE, component="symmetric", weighting="none", norm=2)

aphim <- maor(fitted(rcm), TRUE, component="antisymmetric", weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
aphiu <- maor(fitted(rcu), TRUE, component="antisymmetric", weighting="uniform", norm=2)
aphin <- maor(fitted(rcn), TRUE, component="antisymmetric", weighting="none", norm=2)

cphim <- maor(fitted(rcm), TRUE, TRUE, weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
cphiu <- maor(fitted(rcu), TRUE, TRUE, weighting="uniform", norm=2)
cphin <- maor(fitted(rcn), TRUE, TRUE, weighting="none", norm=2)

maorm <- maor(fitted(rcm), weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
maoru <- maor(fitted(rcu), weighting="uniform", norm=2)
maorn <- maor(fitted(rcn), weighting="none", norm=2)

smaorm <- maor(fitted(rcm), component="symmetric", weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
smaoru <- maor(fitted(rcu), component="symmetric", weighting="uniform", norm=2)
smaorn <- maor(fitted(rcn), component="symmetric", weighting="none", norm=2)

amaorm <- maor(fitted(rcm), component="antisymmetric", weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
amaoru <- maor(fitted(rcu), component="antisymmetric", weighting="uniform", norm=2)
amaorn <- maor(fitted(rcn), component="antisymmetric", weighting="none", norm=2)

cmaorm <- maor(fitted(rcm), cell=TRUE, weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
cmaoru <- maor(fitted(rcu), cell=TRUE, weighting="uniform", norm=2)
cmaorn <- maor(fitted(rcn), cell=TRUE, weighting="none", norm=2)

stopifnot(all.equal(phim, sphim))
stopifnot(all.equal(phiu, sphiu))
stopifnot(all.equal(phin, sphin))
stopifnot(all.equal(phim, sqrt(sum((rcm$assoc$phi)^2))))
stopifnot(all.equal(phiu, sqrt(sum(abs(rcu$assoc$phi)^2))))
stopifnot(all.equal(phin, sqrt(sum(abs(rcn$assoc$phi)^2))))

stopifnot(all(c(aphim, aphiu, aphin) < 1e-10))

stopifnot(all.equal(phim, sqrt(sum(cphim))))
stopifnot(all.equal(phiu, sqrt(sum(cphiu))))
stopifnot(all.equal(phin, sqrt(sum(cphin))))

stopifnot(all(c(amaorm, amaoru, amaorn) - 1 < 1e-10))

stopifnot(all.equal(maorm, exp(sqrt(sum(cmaorm)))))
stopifnot(all.equal(maoru, exp(sqrt(sum(cmaoru)))))
stopifnot(all.equal(maorn, exp(sqrt(sum(cmaorn)))))

stopifnot(all.equal(maorm, smaorm))
stopifnot(all.equal(maoru, smaoru))
stopifnot(all.equal(maorn, smaorn))

# Test on perfectly anti-symmetric association
hmm <- hmskew(ocg1973, nd.symm=0, weighting="marginal", start=NA)
hmu <- hmskew(ocg1973, nd.symm=0, weighting="uniform", start=NA)
hmn <- hmskew(ocg1973, nd.symm=0, weighting="none", start=NA)

stopifnot(all.equal(fitted(hmm), fitted(hmu)))
stopifnot(all.equal(fitted(hmm), fitted(hmn)))

w <- (hmm$assoc$row.weights + hmm$assoc$col.weights)[,1]/2
phim <- maor(fitted(hmm), TRUE, weighting="marginal", norm=2, row.weights=w, col.weights=w)
phiu <- maor(fitted(hmu), TRUE, weighting="uniform", norm=2)
phin <- maor(fitted(hmn), TRUE, weighting="none", norm=2)

sphim <- maor(fitted(hmm), TRUE, component="symmetric", weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
sphiu <- maor(fitted(hmu), TRUE, component="symmetric", weighting="uniform", norm=2)
sphin <- maor(fitted(hmn), TRUE, component="symmetric", weighting="none", norm=2)

aphim <- maor(fitted(hmm), TRUE, component="antisymmetric", weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
aphiu <- maor(fitted(hmu), TRUE, component="antisymmetric", weighting="uniform", norm=2)
aphin <- maor(fitted(hmn), TRUE, component="antisymmetric", weighting="none", norm=2)

cphim <- maor(fitted(hmm), TRUE, TRUE, weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
cphiu <- maor(fitted(hmu), TRUE, TRUE, weighting="uniform", norm=2)
cphin <- maor(fitted(hmn), TRUE, TRUE, weighting="none", norm=2)

maorm <- maor(fitted(hmm), weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
maoru <- maor(fitted(hmu), weighting="uniform", norm=2)
maorn <- maor(fitted(hmn), weighting="none", norm=2)

smaorm <- maor(fitted(hmm), component="symmetric", weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
smaoru <- maor(fitted(hmu), component="symmetric", weighting="uniform", norm=2)
smaorn <- maor(fitted(hmn), component="symmetric", weighting="none", norm=2)

amaorm <- maor(fitted(hmm), component="antisymmetric", weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
amaoru <- maor(fitted(hmu), component="antisymmetric", weighting="uniform", norm=2)
amaorn <- maor(fitted(hmn), component="antisymmetric", weighting="none", norm=2)

cmaorm <- maor(fitted(hmm), cell=TRUE, weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
cmaoru <- maor(fitted(hmu), cell=TRUE, weighting="uniform", norm=2)
cmaorn <- maor(fitted(hmn), cell=TRUE, weighting="none", norm=2)

stopifnot(all.equal(phim, aphim))
stopifnot(all.equal(phiu, aphiu))
stopifnot(all.equal(phin, aphin))
stopifnot(all.equal(phim, sqrt(sum((hmm$assoc$phi)^2))))
stopifnot(all.equal(phiu, sqrt(sum(abs(hmu$assoc$phi)^2))))
stopifnot(all.equal(phin, sqrt(sum(abs(hmn$assoc$phi)^2))))

stopifnot(all(c(sphim, sphiu, sphin) < 1e-10))

stopifnot(all.equal(phim, sqrt(sum(cphim))))
stopifnot(all.equal(phiu, sqrt(sum(cphiu))))
stopifnot(all.equal(phin, sqrt(sum(cphin))))

stopifnot(all(c(smaorm, smaoru, smaorn) - 1 < 1e-10))

stopifnot(all.equal(maorm, exp(sqrt(sum(cmaorm)))))
stopifnot(all.equal(maoru, exp(sqrt(sum(cmaoru)))))
stopifnot(all.equal(maorn, exp(sqrt(sum(cmaorn)))))

stopifnot(all.equal(maorm, amaorm))
stopifnot(all.equal(maoru, amaoru))
stopifnot(all.equal(maorn, amaorn))

# Test on symmetric and anti-symmetric association
# Without starting values, model too often converges to wrong solution
start <- c(6.540,  0.106,  0.407, 0.666, 1.006, -0.581, -0.261, 0.060,
          -4.411, -0.567, -0.310, 0.264, 0.652, -1.794, -1.610, -1.627,
          -0.743, -0.012,  6.311, 0.295, 0.198, -0.015, -0.167, 0.010)
hmm <- hmskew(ocg1973, nd.symm=1, weighting="marginal", start=start)
hmu <- hmskew(ocg1973, nd.symm=1, weighting="uniform", start=start)
hmn <- hmskew(ocg1973, nd.symm=1, weighting="none", start=start)

stopifnot(all.equal(fitted(hmm), fitted(hmu)))
stopifnot(all.equal(fitted(hmm), fitted(hmn)))

w <- (hmm$assoc$row.weights + hmm$assoc$col.weights)[,1]/2
phim <- maor(fitted(hmm), TRUE, weighting="marginal", norm=2, row.weights=w, col.weights=w)
phiu <- maor(fitted(hmu), TRUE, weighting="uniform", norm=2)
phin <- maor(fitted(hmn), TRUE, weighting="none", norm=2)

sphim <- maor(fitted(hmm), TRUE, component="symmetric", weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
sphiu <- maor(fitted(hmu), TRUE, component="symmetric", weighting="uniform", norm=2)
sphin <- maor(fitted(hmn), TRUE, component="symmetric", weighting="none", norm=2)

aphim <- maor(fitted(hmm), TRUE, component="antisymmetric", weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
aphiu <- maor(fitted(hmu), TRUE, component="antisymmetric", weighting="uniform", norm=2)
aphin <- maor(fitted(hmn), TRUE, component="antisymmetric", weighting="none", norm=2)

cphim <- maor(fitted(hmm), TRUE, TRUE, weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
cphiu <- maor(fitted(hmu), TRUE, TRUE, weighting="uniform", norm=2)
cphin <- maor(fitted(hmn), TRUE, TRUE, weighting="none", norm=2)

maorm <- maor(fitted(hmm), weighting="marginal", norm=2,
              row.weights=w, col.weights=w)
maoru <- maor(fitted(hmu), weighting="uniform", norm=2)
maorn <- maor(fitted(hmn), weighting="none", norm=2)

smaorm <- maor(fitted(hmm), component="symmetric", weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
smaoru <- maor(fitted(hmu), component="symmetric", weighting="uniform", norm=2)
smaorn <- maor(fitted(hmn), component="symmetric", weighting="none", norm=2)

amaorm <- maor(fitted(hmm), component="antisymmetric", weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
amaoru <- maor(fitted(hmu), component="antisymmetric", weighting="uniform", norm=2)
amaorn <- maor(fitted(hmn), component="antisymmetric", weighting="none", norm=2)

cmaorm <- maor(fitted(hmm), cell=TRUE, weighting="marginal", norm=2,
               row.weights=w, col.weights=w)
cmaoru <- maor(fitted(hmu), cell=TRUE, weighting="uniform", norm=2)
cmaorn <- maor(fitted(hmn), cell=TRUE, weighting="none", norm=2)

stopifnot(all.equal(phim, sqrt(sphim^2 + aphim^2)))
stopifnot(all.equal(phiu, sqrt(sphiu^2 + aphiu^2)))
stopifnot(all.equal(phin, sqrt(sphin^2 + aphin^2)))
stopifnot(all.equal(sphim, sqrt(sum((hmm$assoc$phi)^2))))
stopifnot(all.equal(sphiu, sqrt(sum(abs(hmu$assoc$phi)^2))))
stopifnot(all.equal(sphin, sqrt(sum(abs(hmn$assoc$phi)^2))))
stopifnot(all.equal(aphim, sqrt(sum((hmm$assoc.hmskew$phi)^2))))
stopifnot(all.equal(aphiu, sqrt(sum(abs(hmu$assoc.hmskew$phi)^2))))
stopifnot(all.equal(aphin, sqrt(sum(abs(hmn$assoc.hmskew$phi)^2))))

stopifnot(all.equal(phim, sqrt(sum(cphim))))
stopifnot(all.equal(phiu, sqrt(sum(cphiu))))
stopifnot(all.equal(phin, sqrt(sum(cphin))))

stopifnot(all.equal(maorm, exp(sqrt(sum(cmaorm)))))
stopifnot(all.equal(maoru, exp(sqrt(sum(cmaoru)))))
stopifnot(all.equal(maorn, exp(sqrt(sum(cmaorn)))))


# Test for phi computed from UNIDIFF two-way interaction coefficients
data(yaish)
tab <- aperm(yaish[,,-7], 3:1)

# 1-norm
# Currently disabled until supported correctly
# u1m <- unidiff(tab, weighting="marginal", norm=1)
# rp <- prop.table(margin.table(tab, 1))
# cp <- prop.table(margin.table(tab, 2))
# stopifnot(all.equal(u1m$unidiff$phi,
#                     maor(fitted(u1m)[,,1], TRUE, weighting="marginal", norm=1, rp, cp)))
# stopifnot(all.equal(u1m$unidiff$phi,
#                     maor(fitted(u1m)[,,1], TRUE, weighting="marginal", norm=1, rp, cp)))
# 
# u1u <- unidiff(tab, weighting="uniform", norm=1)
# stopifnot(all.equal(u1u$unidiff$phi,
#                     maor(fitted(u1u)[,,1], TRUE, weighting="uniform", norm=1)))
# stopifnot(all.equal(u1u$unidiff$phi * exp(u1u$unidiff$layer$qvframe$estimate[2]),
#                     maor(fitted(u1u)[,,2], TRUE, weighting="uniform", norm=1)))
# 
# u1n <- unidiff(tab, weighting="none", norm=1)
# stopifnot(all.equal(u1n$unidiff$phi,
#                     maor(fitted(u1n)[,,1], TRUE, weighting="none", norm=1)))
# stopifnot(all.equal(u1n$unidiff$phi * exp(u1u$unidiff$layer$qvframe$estimate[2]),
#                     maor(fitted(u1n)[,,2], TRUE, weighting="none", norm=1)))


# 2-norm
u2m <- unidiff(tab, weighting="marginal", norm=2)
stopifnot(all.equal(u2m$unidiff$phi,
                    maor(fitted(u2m), TRUE, norm=2)[1], weighting="marginal", check.attributes=FALSE))
stopifnot(all.equal(u2m$unidiff$phi * exp(u2m$unidiff$layer$qvframe$estimate[2]),
                    maor(fitted(u2m), TRUE, norm=2)[2], weighting="marginal", check.attributes=FALSE))

u2u <- unidiff(tab, weighting="uniform", norm=2)
stopifnot(all.equal(u2u$unidiff$phi,
                    maor(fitted(u2u), TRUE, weighting="uniform", norm=2)[1], check.attributes=FALSE))
stopifnot(all.equal(u2u$unidiff$phi * exp(u2u$unidiff$layer$qvframe$estimate[2]),
                    maor(fitted(u2u), TRUE, weighting="uniform", norm=2)[2], check.attributes=FALSE))

u2n <- unidiff(tab, weighting="none", norm=2)
stopifnot(all.equal(u2n$unidiff$phi,
                    maor(fitted(u2n), TRUE, weighting="none", norm=2)[1], check.attributes=FALSE))
stopifnot(all.equal(u2n$unidiff$phi * exp(u2n$unidiff$layer$qvframe$estimate[2]),
                    maor(fitted(u2n), TRUE, weighting="none", norm=2)[2], check.attributes=FALSE))

stopifnot(all.equal(maor(yaish),
                    apply(yaish, 3, maor,
                          row.weights=margin.table(yaish, 1),
                          col.weights=margin.table(yaish, 2))))


###
## Comparison with mean/sum of all spanning odds ratios
###
# Can be set to arbitrary values
nr <- 4
nc <- 5

norm <- 2

or <- function(tab) {
    or1 <- function(i, j, i2, j2) (tab[i, j] * tab[i2, j2]) / (tab[i, j2] * tab[i2, j])
    or <- array(NA, c(nrow(tab), ncol(tab), nrow(tab), ncol(tab)))

    for(i in 1:nrow(tab))
      for(j in 1:ncol(tab))
        for(i2 in 1:nrow(tab))
          for(j2 in 1:ncol(tab))
            if(i2 != i && j2 != j)
              or[i, j, i2, j2] <- or1(i, j, i2, j2)
    or
}

wlor2 <- function(tab) {
    rp <- prop.table(margin.table(tab, 1)) * nrow(tab)
    cp <- prop.table(margin.table(tab, 2)) * ncol(tab)
    wlor2 <- w <- array(NA, c(nrow(tab), ncol(tab), nrow(tab), ncol(tab)))

    for(i in 1:nrow(tab))
      for(j in 1:ncol(tab))
        for(i2 in 1:nrow(tab))
          for(j2 in 1:ncol(tab))
            if(i2 != i && j2 != j) {
              wlor2[i, j, i2, j2] <- log((tab[i, j] * tab[i2, j2]) / (tab[i, j2] * tab[i2, j]))^norm
              w[i, j, i2, j2] <- rp[i] * cp[j] * rp[i2] * cp[j2]
            }

    wlor2 * w / sum(w, na.rm=TRUE)
}

## Unweighted mean
# General case
res <- replicate(10, {
    tab <- matrix(rpois(nr*nc, 1000), nr, nc) + .5
    rp <- rep(1/nr, nr)
    cp <- rep(1/nc, nc)
    c(maor(tab, weighting="uniform", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(mean(abs(log(or(tab)))^norm, na.rm=TRUE)^(1/norm)))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))

# 2x2 table (equality with single odds ratio)
res <- replicate(10, {
    tab <- matrix(rpois(2*2, 1000), 2, 2) + .5
    rp <- cp <- rep(1/2, 2)
    c(maor(tab, weighting="uniform", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(sqrt(log(tab[1,1] * tab[2,2] / (tab[1,2] * tab[2,1]))^2)))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))


## Unweighted sum
# General case
res <- replicate(10, {
    tab <- matrix(rpois(nr*nc, 1000), nr, nc) + .5
    rp <- rep(1, nr)
    cp <- rep(1, nc)
    c(maor(tab, weighting="none", norm=norm),
      exp(sum(abs(log(or(tab)))^norm, na.rm=TRUE)^(1/norm)))
})

stopifnot(all.equal(res[1,], res[2,]))

# 2x2 table (equality with single odds ratio)
res <- replicate(10, {
    tab <- matrix(rpois(2*2, 1000), 2, 2) + .5
    rp <- cp <- rep(1, 2)
    c(maor(tab, weighting="none", norm=norm),
      exp(sqrt(4 * log(tab[1,1] * tab[2,2] / (tab[1,2] * tab[2,1]))^2)))
})

stopifnot(all.equal(res[1,], res[2,]))


## Marginal-weighted mean
# General case
res <- replicate(10, {
    tab <- matrix(rpois(nr*nc, 100), nr, nc) +.5
    rp <- prop.table(margin.table(tab, 1))
    cp <- prop.table(margin.table(tab, 2))
    c(maor(tab, weighting="marginal", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(abs(sum(wlor2(tab), na.rm=TRUE))^(1/norm)))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))

# 2x2 table (equality with single odds ratio)
res <- replicate(10, {
    tab <- matrix(rpois(2*2, 100), 2, 2) + .5
    rp <- prop.table(margin.table(tab, 1))
    cp <- prop.table(margin.table(tab, 2))
    c(maor(tab, weighting="marginal", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(abs(log(tab[1,1] * tab[2,2] / (tab[1,2] * tab[2,1])))))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))
