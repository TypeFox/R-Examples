
# mlplot demo
# rv package
# 2007-07-07
#

#
# prepare a fictional data set
#
n.rows <- 4; n.cols <- 5; n <- (n.rows*n.cols)
# Draw some fixed numbers
mu.true <- rnorm(1:n.rows, mean=1:n.rows, sd=1)
sigma.true <- 1.0
theta <- rvmatrix(rvnorm(n=n.cols, mean=mu.true, sd=sigma.true), nrow=n.rows)
rvcols <- c("red", "blue", "green", "orange")
#
col.labels <- paste("Time", 1:n.cols, sep=":")
row.labels <- paste("Unit", 1:n.rows, sep=":")
dimnames(theta) <- list(row.labels, col.labels)
#
# mlplot examples
#
par(mfrow=c(2,2))
mlplot(theta, 
  main=paste(n.cols,"repeated measurements for each of the", n.rows,
    "units\ncyan: row means; vertical line: E(mean(theta))")
)
#
grand.mean.estimate <- E(mean(theta))
abline(v=grand.mean.estimate)
abline(v=0, lty="dotted")
theta.row.means <- apply.rv(theta, 1, mean)
mlplot(theta.row.means, pch=19, y.shift=(n.cols+1)/20, add=TRUE, rvcol="cyan")
#
unit.colors <- paste(1:n.rows, rvcols, sep=":", collapse=",")
mlplot(theta, main=paste("Each unit is colored", unit.colors, sep="\n"), rvcol=rvcols) 
abline(v=0, lty="dotted")
mlplot(theta.row.means, pch=19, y.shift=(n.cols+1)/20, add=TRUE, rvcol="cyan")
#
mlplot(t(theta), main="Same random matrix, transposed\npurple: column means")
abline(v=0, lty="dotted")
abline(v=grand.mean.estimate)
theta.col.means <- apply.rv(t(theta), 1, mean)
mlplot(theta.col.means, pch=19, y.shift=(n.cols+1)/20, add=TRUE, rvcol="purple")
#
Cols <- matrix(rvcols, nrow=n.rows, ncol=n.cols)
mlplot(t(theta), main=paste("Transposed matrix, each unit is colored", unit.colors, sep="\n"), rvcol=t(Cols))
abline(v=0, lty="dotted")
abline(v=grand.mean.estimate)
mlplot(theta.col.means, pch=19, y.shift=(n.cols+1)/20, add=TRUE, rvcol="purple")


# end mlplotdemo.R

