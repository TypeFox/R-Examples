## ------------------------------------------------------------------------
set.seed(1)

p     <- 500
probs <- runif(p, 0.1, 0.5)
probs <- t(probs) %x% matrix(1,p,2)
X0    <- matrix(rbinom(2*p*p, 1, probs), p, 2*p)
X     <- X0 %*% (diag(p) %x% matrix(1,2,1))

## ---- results = "asis"---------------------------------------------------
pander::pandoc.table(X[1:10, 1:10])

## ------------------------------------------------------------------------
group <- c(rep(1:20, each=3),
           rep(21:40, each=4),
           rep(41:60, each=5),
           rep(61:80, each=6),
           rep(81:100, each=7))
group.id <- grpSLOPE::getGroupID(group)
n.group <- length(group.id)
group.length <- sapply(group.id, FUN=length)

## ------------------------------------------------------------------------
ind.relevant <- sample(1:n.group, 10)
print(sort(ind.relevant))

## ------------------------------------------------------------------------
b <- rep(0, p)
for (j in ind.relevant) {
  # generate effect sizes from the Uniform(0,1) distribution
  b[group.id[[j]]] <- runif(group.length[j])
}

# generate the response vector
y <- X %*% b + rnorm(p)

## ------------------------------------------------------------------------
library(grpSLOPE)

result <- grpSLOPE::grpSLOPE(X=X, y=y, group=group, fdr=0.1)

## ------------------------------------------------------------------------
result$selected

## ------------------------------------------------------------------------
# estimated sigma (true sigma is equal to one)
result$sigma
# first 14 entries of b estimate:
result$beta[1:14]

## ------------------------------------------------------------------------
plot(result$lambda[1:10], xlab = "Index", ylab = "Lambda", type="l")

## ---- results = "asis"---------------------------------------------------
n.selected    <- length(result$selected)
true.relevant <- names(group.id)[ind.relevant]
truepos       <- intersect(result$selected, true.relevant)

n.truepos  <- length(truepos)
n.falsepos <- n.selected - n.truepos

gFDP <- n.falsepos / max(1, n.selected)
pow <- n.truepos / length(true.relevant)

print(paste("gFDP =", gFDP))
print(paste("Power =", pow))

