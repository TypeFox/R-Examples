mantel <- function(formula = formula(data), data = sys.parent(), nperm = 1000, mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95)
{
# Mantel test 
# Written by Sarah C. Goslee
# 27 June 2000
# Updated 5 April 2001
#
# formula is y ~ x + n1 + n2 + n3 + ...
# NOT y ~ x | n1 + n2 + n3 + ... 
# The | means something else in S-Plus formulas.
#
# Uses C for permutation and bootstrap routines.
#
# This version calculates partial coefficients by permuting the y matrix.
#
# Will calculate the simple correlation or n-th order partial correlation  
# between two distance matrices in either of two ways: Pearson (mrank=FALSE) 
# or Spearman (mrank=TRUE)
#
# A permutation test is used to calculate the significance of r.
# The permutation test was designed to be relatively fast, but because of the
# way this was done, there is a possibility of repeating permutations of
# 1/n! where the distance matrix is n by n. In particular, for small matrices 
# n < 8 or so, it may be better to enumerate the permutations. 
#
#
# As an added bonus, this function offers the option of calculating 
# bootstrapped confidence limits for the correlation coefficient.
# nboot is the number of iterations.
# pboot is the level to resample at.
# cboot is the desired confidence limit.
# NOTE: This is not bootstrapping with replacement. That doesn't make
# much sense for dissimilarities because of the possibility of duplicates.
# The dissimilarity between a sample and itself is always zero.
#
# mantel returns a five-element list:
# mantelr is the correlation.
# pval1 is the one-sided p-value (null hypothesis r <= 0) (0 if nperm == 0).
# pval2 is the one-sided p-value (null hypothesis r >= 0) (0 if nperm == 0).
# pval3 is the two-sided p-value (null hypothesis r = 0) (0 if nperm == 0).
# llim is the lower confidence limit.
# ulim is the upper confidence limit.
# 
# requires mantel.c (Included in ecodist.c.)
#

# Stuff R needs to be able to use a formula
        m <- match.call(expand.dots = FALSE)
        m2 <- match(c("formula", "data"), names(m), nomatch=0)
        m <- m[c(1, m2)]
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        m <- as.matrix(m)

# End of R stuff. m is now the data for the Mantel test as
# columns y, x, n1, n2, n3, ...
# Determine the size of the matrices & do some error checking.
        n <- (1 + sqrt(1 + 8 * nrow(m)))/2
        if(abs(n - round(n)) > 0.0000001)
		stop("Matrix not square.\n")
	n <- round(n)
        if(ncol(m) < 2) stop("Not enough data. \n")

# If there are only x and y, then use the data as is.
        if(dim(m)[[2]] == 2) {
                ymat <- as.vector(m[, 1])
                xmat <- as.vector(m[, 2])
                if(mrank) {
                        ymat <- rank(ymat)
                        xmat <- rank(xmat)
                }
                ycor <- ymat
                xcor <- xmat
        }
        else {
                ymat <- as.vector(m[, 1])
                omat <- m[, -1]
                if(mrank) {
                        ymat <- rank(ymat)
                        omat <- apply(omat, 2, rank)
                }
                omat <- cbind(rep(1, length(ymat)), omat)
                xmat <- as.vector(omat[, 2])
                omat <- omat[, -2]
		omat <- as.matrix(omat)
                ycor <- lm.fit(omat, ymat)$residuals
                xcor <- lm.fit(omat, xmat)$residuals
        }

        mantelr <- cor(xcor, ycor)

# Convert matrices to column order for compatibility with C routines.

        xmat <- full(xmat)
        ymat <- full(ymat)
        xmat <- xmat[col(xmat) > row(xmat)]
        ymat <- ymat[col(ymat) > row(ymat)]

        if(dim(m)[[2]] > 2) {
                for(i in 2:dim(omat)[[2]]) {
                        curcoll <- omat[, i]
                        curcoll <- full(curcoll)
                        curcoll <- curcoll[col(curcoll) > row(curcoll)]
                        omat[, i] <- curcoll
                }
        }
# If using a permutation test, start here:

        if(nperm > 0) {

# Set up the arrays needed.
                zstats <- numeric(nperm)
                tmat <- matrix(0, n, n)
                rarray <- rep(0, n)
                if(dim(m)[[2]] == 2) {

# Standardize the columns of the matrices so
# that z = r and we can do 2-tailed tests.
 
			ncor <- length(xmat)

			w1 <- sum(xmat)/ncor
			w2 <- sum(xmat^2)
			w2 <- sqrt(w2/ncor - w1^2)
			xmat <- (xmat - w1)/w2

			w1 <- sum(ymat)/ncor
			w2 <- sum(ymat^2)
			w2 <- sqrt(w2/ncor - w1^2)
			ymat <- (ymat - w1)/w2

             cresults <- .C("permute",
                                as.double(xmat),
                                as.double(ymat),
                                as.integer(n),
                                as.integer(length(xmat)),
                                as.integer(nperm),
                                zstats = as.double(zstats),
                                as.double(as.vector(tmat)),
                                as.integer(rarray),
										  PACKAGE = "ecodist")
		}
                else {
			tomat <- t(omat)
                        hmat <- solve(tomat %*% omat)
                        hmat <- hmat %*% tomat
			bmat <- rep(0, ncol(omat))

                        xcor <- as.vector(lm.fit(omat, xmat)$residuals)
                        ycor <- as.vector(lm.fit(omat, ymat)$residuals)

# Standardize the columns of the matrices so
# that z = r and we can do 2-tailed tests.
 
			ncor <- length(xcor)

			w1 <- sum(xcor)/ncor
			w2 <- sum(xcor^2)
			w2 <- sqrt(w2/ncor - w1^2)
			xcor <- (xcor - w1)/w2

			w1 <- sum(ycor)/ncor
			w2 <- sum(ycor^2)
			w2 <- sqrt(w2/ncor - w1^2)
			ycor <- (ycor - w1)/w2

         cresults <- .C("permpart",
                                as.double(as.vector(hmat)),
                                bmat = as.double(bmat),
				as.double(as.vector(omat)),
                                as.double(ymat),
                                as.double(xcor),
                                ycor = as.double(ycor),
                                as.integer(n),
				as.integer(length(bmat)),
                                as.integer(length(xmat)),
                                as.integer(nperm),
                                zstats = as.double(zstats),
                                as.double(as.vector(tmat)),
                                as.integer(rarray),
										  PACKAGE = "ecodist")
                }
                zstats <- cresults$zstats       # Calculate the p-values.
                pval1 <- length(zstats[zstats >= zstats[1]])/nperm
                pval2 <- length(zstats[zstats <= zstats[1]])/nperm
                pval3 <- length(zstats[abs(zstats) >= abs(zstats[1])])/nperm
        }
        else {
                pval1 <- 0
                pval2 <- 0
                pval3 <- 0
        }
# If not using a permutation test, return 0 for the p-values.
        if(nboot > 0) {
                if(dim(m)[[2]] == 2) {
                        ycor <- ymat
                        xcor <- xmat
                }
                else {
                        xcor <- as.vector(lm.fit(omat, xmat)$residuals)
                        ycor <- as.vector(lm.fit(omat, ymat)$residuals)
                }
                bootcor <- numeric(nboot)
                rarray <- numeric(n)
                rmat <- numeric(length(xcor))
                xdif <- numeric(length(xcor))
                ydif <- numeric(length(xcor))
                cresults <- .C("bootstrap",
                        as.double(xcor),
                        as.double(ycor),
                        as.integer(n),
                        as.integer(length(xcor)),
                        as.integer(nboot),
                        as.double(pboot),
                        bootcor = as.double(bootcor),
                        as.integer(rarray),
                        as.integer(rmat),
                        as.double(xdif),
                        as.double(ydif),
								PACKAGE = "ecodist")
                bootcor <- cresults$bootcor
                bootcor <- sort(bootcor)
                pval <- (1 - cboot)/2
                llim <- quantile(bootcor, pval)
                ulim <- quantile(bootcor, 1 - pval)
        }
        else {
                llim <- 0
                ulim <- 0
        }
        c(mantelr = mantelr, pval1 = pval1, pval2 = pval2, pval3 = pval3, llim = llim, ulim = ulim)
}
