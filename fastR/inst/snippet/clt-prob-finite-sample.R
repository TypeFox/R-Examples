x <- c(1,2,4,4,9)
mu <- sum(x * 0.2); mu                 # population mean
v <- sum(x^2 *0.2) - mu^2; v           # population variance
pairsums <- outer(x,x,"+")             # compute 25 sums
pairmeans <- pairsums/2

# sampling distribution with SRS
srs.means <- as.vector(pairmeans[lower.tri(pairmeans)]); srs.means
iid.means <- as.vector(pairmeans); iid.means

srs.mean <- sum(srs.means * 0.1); srs.mean
srs.var <- sum(srs.means^2 * 0.1) - srs.mean^2; srs.var
v/2 * (5-2) / (5-1)
sqrt(v/2 * (5-2) / (5-1))

var(srs.means)   # INCORRECT variance

# sampling distribution with iid sample
iid.mean <- sum(iid.means * 0.04); iid.mean
iid.var <- sum(iid.means^2 * 0.04) - iid.mean^2; iid.var
v/2
sqrt(v/2)

var(iid.means)   # INCORRECT variance
