
##############################################################

# A high-degree polynomial system (R.B. Kearfoot, ACM 1987)
# There are 12 real roots (and 126 complex roots to this system!)

##############################################################

require("BB")

hdp <- function(x) {
f <- rep(NA, length(x))
f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
f[3] <- x[1]^2 + x[2]^2 - 0.265625
f
}

# Multiple starting values
set.seed(123)
p0 <- matrix(runif(600), 200, 3)  # 200 starting values, each of length 3

# dfsane() with default parameters
ans.df <- matrix(NA, nrow(p0), ncol(p0))
for (i in 1:nrow(p0)) {
tmp <- dfsane(par=p0[i, ], fn=hdp, control=list(trace=FALSE))
if (tmp$conv == 0) ans.df[i, ] <- tmp$par	
}

ans.df <- ans.df[!is.na(ans.df[,1]), ]
dim(ans.df) # note that only 69 successes

# BBsolve() 
ans <- multiStart(par=p0, fn=hdp)
sum(ans$conv)  # number of successful runs = 190
pmat <- ans$par[ans$conv, ] # selecting only converged solutions
pc <- princomp(pmat)

biplot(pc)  # you can see all 12 solutions beautifully like on a clock!

