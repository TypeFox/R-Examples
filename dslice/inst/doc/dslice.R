### R code from vignette source 'dslice.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(dslice)


###################################################
### code chunk number 2: generate_data
###################################################
n <- 100
mu <- 0.5
y <- c(rnorm(n, -mu, 1), rnorm(n, mu, 1))
x <- c(rep("G1", n), rep("G2", n))


###################################################
### code chunk number 3: convert
###################################################
x <- relabel(x)


###################################################
### code chunk number 4: get_dim
###################################################
xdim <- max(x) + 1


###################################################
### code chunk number 5: order_x
###################################################
x <- x[order(y)]


###################################################
### code chunk number 6: run_ds_k
###################################################
lambda <- 1.0
dsres <- ds_k(x, xdim, lambda, slice = TRUE)
dsres


###################################################
### code chunk number 7: illustrate_slice
###################################################
colnames(dsres$slices) <- c("G1", "G2", "total")
slice_show(dsres$slices)


###################################################
### code chunk number 8: use_ds_eqp_k
###################################################
n <- 100
mu <- 0.5
y <- c(rnorm(n, -mu, 1), rnorm(n, mu, 1))
x <- c(rep("1", n), rep("2", n))
x <- relabel(x)
x <- x[order(y)]
xdim <- max(x) + 1
lambda <- 1.0
dsres <- ds_eqp_k(x, xdim, lambda, slice = TRUE)  


###################################################
### code chunk number 9: data_for_ds_test
###################################################
n <- 100
y <- c(rnorm(n, -mu, 1), rnorm(n, mu, 1))

##  generate x in this way:
x <- c(rep(0, n), rep(1, n))
x <- as.integer(x)

##  or in this way:
x <- c(rep("G1", n), rep("G2", n))
x <- relabel(x)


###################################################
### code chunk number 10: ds_test
###################################################
lambda <- 1.0
dsres <- ds_test(y, x, type = "eqp", lambda = 1, rounds = 100)


###################################################
### code chunk number 11: one_sample_test
###################################################
##  One-sample test
n <- 100
mu <- 0.5
y <- rnorm(n, mu, 1)
lambda <- 1.0
alpha <- 1.0
dsres <- ds_test(y, "pnorm", 0, 1, lambda = 1, alpha = 1, rounds = 100)
dsres <- ds_test(y, "pnorm", 0, 1, type = "ds", lambda = 1, alpha = 1)
dsres <- ds_test(y, "pnorm", 0, 1, type = "eqp", lambda = 1, rounds = 100)
dsres <- ds_test(y, "pnorm", 0, 1, type = "eqp", lambda = 1)


###################################################
### code chunk number 12: ds_1
###################################################
n <- 100
mu <- 0.5
x <- rnorm(n, mu, 1)
y <- pnorm(sort(x), 0, 1) 
lambda <- 1.0
alpha <- 1.0
dsres <- ds_1(y, lambda, alpha)
dsres <- ds_eqp_1(y, lambda)


###################################################
### code chunk number 13: see_data
###################################################
data()


###################################################
### code chunk number 14: load_gsa
###################################################
data(gsa_exp)
data(gsa_label)
data(gsa_set)


###################################################
### code chunk number 15: rank_function
###################################################
fc <- function(mat, label)
{
  d0 <- apply(x[,which(label == 0)], 1, mean)
  d1 <- apply(x[,which(label == 1)], 1, mean)
  d <- d1 / d0
  return(order(d))
}


###################################################
### code chunk number 16: diffslices
###################################################
set.seed(1)
n<-30
x1<-rep(1,n)
x2<-rep(2,n)
y1<-c(rnorm(n/2,-2,0.5),rnorm(n/2,2,0.5))
y2<-rnorm(n,0,1)

xrange <- c(0.5,2.5)
yrange <- c(-3.2,3.2)

par(mfrow=c(1,2))

plot(x1, y1, xlim=xrange, ylim=yrange, pch=paste(1), col=2, xlab="X", ylab="Y", xaxt="n", cex=1.5)
axis(side=1, at=c(1,2), label=c(1,2))
par(new=T)
plot(x2, y2, xlim=xrange, ylim=yrange, pch=paste(2), col=4, xlab="", ylab="", xaxt="n", cex=1.5)
axis(side=1, at=c(1,2), label=c(1,2))
abline(h=0, lty=2, lwd=2)

plot(x1, y1, xlim=xrange, ylim=yrange, pch=paste(1), col=2, xlab="X", ylab="Y", xaxt="n", cex=1.5)
axis(side=1, at=c(1,2), label=c(1,2))
par(new=T)
plot(x2, y2, xlim=xrange, ylim=yrange, pch=paste(2), col=4, xlab="", ylab="", xaxt="n", cex=1.5)
axis(side=1, at=c(1,2), label=c(1,2))
abline(h=1, lty=2, lwd=2)
abline(h=-1, lty=2, lwd=2)


