library(smaa)

# First alternative beats second on all criteria
# --> NA central weights
N <- 1E4; m <- 2; n <- 3
mnames <- c("A", "B")
nnames <- c("a", "b", "c")
meas <- array(dim=c(N,m,n))
meas[,1,] <- 1
meas[,2,] <- 0
dimnames(meas) <- list(NULL, mnames, nnames)
pref <- matrix(1/3, nrow=N, ncol=n)
dimnames(pref) <- list(NULL, nnames)

centroid <- rep(1/3, 3)
names(centroid) <- nnames

result <- smaa(meas, pref)
stopifnot(all.equal(result$cw[1,], centroid))
stopifnot(all(is.na(result$cw[2,])))
stopifnot(all.equal(rownames(result$cw), mnames))

ranks <- smaa.ranks(smaa.values(meas, pref))
cw <- smaa.cw(ranks, pref)
stopifnot(all.equal(cw[1,], centroid))
stopifnot(all(is.na(cw[2,])))
stopifnot(all.equal(rownames(cw), mnames))

cf <- smaa.cf(meas, cw)
stopifnot(all.equal(unname(cf$cf), c(1, NA)))
stopifnot(all.equal(names(cf$cf), mnames))


# First alternative beats second in a single iteration
# --> R gives us a vector instead of a matrix (gee, thanks!)
meas[1,1,] <- 0
meas[1,2,] <- 1

result <- smaa(meas, pref)
stopifnot(all.equal(result$cw[1,], centroid))
stopifnot(all.equal(result$cw[2,], centroid))

ranks <- smaa.ranks(smaa.values(meas, pref))
cw <- smaa.cw(ranks, pref)
stopifnot(all.equal(cw[1,], centroid))
stopifnot(all.equal(cw[2,], centroid))
