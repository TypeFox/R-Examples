init.rand <- function(x,y,K,intercept) {
#
# Function init.rand.  To form starting values for the EM
# algorithm, randomly.
#

tmp <- lm(y~x-1)
ccc <- coef(tmp)
ncc <- length(ccc)
sdc <- 0.3*(abs(ccc))
vvv <- summary(tmp)$sigma**2
lll <- 1/K

rslt <- list()
for(j in 1:K) {
	cft <- rnorm(ncc,ccc,sdc)
	rslt[[j]] <- list(beta=cft,sigsq=vvv,lambda=lll)
}
rslt
}
