
library("mlt")
library("np")
set.seed(29)

### true dgp
rY <- function(n) rgamma(n, shape = 5)
pY <- function(x) pgamma(x, shape = 5)
dY <- function(x) dgamma(x, shape = 5)

### generate y and set-up basis
y <- sort(rY(100))
Bb <- Bernstein_basis(numeric_var("y", support = c(0, max(y) + .1)), order = 10,
                      ui = "increasing")

mydata <- data.frame(y = y)
opt <- mlt(ctm(response = Bb), data = mydata)
d <- opt$todistr

### evaluate on grid
yn <- mkgrid(Bb, n = 50)$y
### eval estimated h and h'
h <- predict(Bb, newdata = data.frame(y = yn), coef = opt$par)
h1 <- predict(Bb, newdata = data.frame(y = yn), deriv = c(y = 1), coef = opt$par)

### plot
layout(matrix(1:2, ncol = 1))
plot(yn, d$p(h), type = "l", col = "red", main = "distribution")
lines(yn, pY(yn), col = "blue")
lines(ecdf(y), col = "grey", cex = .1)
lines(yn, predict(npudist(npudistbw(~ y)), newdata = data.frame(y = yn)), 
      col = "magenta")
rug(y)
plot(yn, d$d(h) * h1,  type = "l", col = "red", main = "density")
lines(yn, dY(yn), col = "blue")
lines(yn, predict(npudens(npudensbw(~ y)), newdata = data.frame(y = yn)), 
      col = "magenta")
lines(density(y), col = "darkgreen")
legend("topright", lwd = 1, legend = c("trafo", "true", "np", "density"),
       col = c("red", "blue", "magenta", "darkgreen"))

