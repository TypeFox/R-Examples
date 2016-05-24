# Tests using the arsenate data and Deming regression
library(deming)
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

# This verifies that the data set is correct
lfit1 <- lm(aes~aas, arsenate)
lfit2 <- lm(aas~aes, arsenate)
aeq(lfit1$coef, c(.544, .845), tol=.001)  # results from the Ripley paper
aeq(lfit2$coef, c(-.299, 1.089), tol=.001)

wfit1 <- lm(aes~aas, arsenate, weight=1/se.aes^2)
wfit2 <- lm(aas~aes, arsenate, weight=1/se.aas^2)
aeq(wfit1$coef, c(0.005, .890), tol=.001) # results from the Ripley paper
aeq(wfit2$coef, c(-.167, 0.631), tol=.001)

#
# Check the jackknife results directly
#
dfit <- deming(aes ~ aas, data=arsenate, xstd=se.aas, ystd=se.aes, dfbeta=TRUE)
dtest <- matrix(0, nrow=nrow(arsenate), ncol=2)
for (i in 1:nrow(dtest)) {
    tfit <- deming(aes ~ aas, data=arsenate, xstd=se.aas, ystd=se.aes, 
                   subset=(-i))
    dtest[i,] <- coef(dfit) - coef(tfit)
}
# The Deming computation uses optimize, whose default tolerance is
#  .Machine$double.eps^.25, in turn leading to less precision in the 
#  replication due to different starting values.
all.equal(dfit$dfbeta, dtest, tol=1e-5)
all.equal(crossprod(dfit$dfbeta), dfit$var)
aeq(coef(dfit), c(.106448, 0.973000), tol=1e-5) #results in paper


# scaled residuals
wt <- 1/(arsenate$se.aes^2 + (dfit$coef[2] * arsenate$se.aas)^2)
rr <- with(arsenate, (aes- (dfit$coef[1] + dfit$coef[2]*aas))*sqrt(wt))
# There is an error in the paper here; it claims to divide by n-2 =28 but 
#  the printed numeric value divides by n.
aeq(mean(rr^2), 1.27, tol=.01)

wtx <- sum(wt*arsenate$aas)/sum(wt)
seb <- 1/sqrt(sum(wt * (arsenate$aas-wtx)^2))  #claimed se for slope
sea <- with(arsenate, sqrt(sum(wt*aas^2)/ (sum(wt)*sum(wt*(aas-wtx)^2))))
aeq(c(sea, seb), dfit$se)

