library("aroma.core")

# Sample scatter data
n <- 10e3
x <- rnorm(n=n)
y <- rnorm(n=n)
xy <- cbind(x=x, y=sin(x)+y/5)

# Bin data and estimate densities
xyd <- binScatter(xy)

layout(matrix(1:4, nrow=2))
par(mar=c(5,4,2,1))

# Plot data
plot(xyd, pch=1)

# Thin scatter data by subsampling
rhos <- c(1/3, 1/4, 1/6)
for (kk in seq(along=rhos)) {
  xyd2 <- subsample(xyd, size=rhos[kk])
  points(xyd2, pch=1, col=kk+1)
}

for (kk in seq(along=rhos)) {
  xyd2 <- subsample(xyd, size=rhos[kk])
  plot(xyd2, pch=1, col=kk+1)
  mtext(side=3, line=0, sprintf("Density: %.1f%%", 100*rhos[kk]))
}
