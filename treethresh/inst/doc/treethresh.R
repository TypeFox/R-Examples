### R code from vignette source 'treethresh.Rnw'

###################################################
### code chunk number 1: treethresh.Rnw:13-14
###################################################
library(treethresh)


###################################################
### code chunk number 2: treethresh.Rnw:78-92
###################################################
from <- -4
to <- 4
x <- seq(from,to,length.out=300)
w <- 0.3
t <- tfromw(w)
plot(x,x,type="n",xaxt="n", yaxt="n",xlab="x",ylab=expression(mu(x)))
usr <- par()$usr
rect(-t,usr[3]+0.02,t,usr[4]-0.02,col=grey(0.85),border=NA)
lines(x,postmed.laplace(x,w=w))
lines(c(from,-t,NA,-t,t,NA,t,to),c(from,-t,NA,0.02,0.02,NA,t,to),lty=2)
lines(c(from,-t,t,to),c(from+t,-0.02,-0.02,to-t),lty=3)
axis(1,at=c(-t,0,t),labels=c(expression(-t[w]),0,expression(t[w])))
axis(2,at=0)
legend("topleft",lty=1:3,legend=c("post.med.","hard","soft"))


###################################################
### code chunk number 3: treethresh.Rnw:156-157
###################################################
w.true <- c(rep(0.15,400),rep(0.6,300),rep(0.05,300))


###################################################
### code chunk number 4: treethresh.Rnw:162-164
###################################################
# Set the seed so that the vignettes are all the same ...
set.seed(111)


###################################################
### code chunk number 5: treethresh.Rnw:167-173
###################################################
mu <- numeric(length(w.true))
non.zero.entry <- runif(length(mu))<w.true
num.non.zero.entries <- sum(non.zero.entry)
mu[non.zero.entry] <- rexp(num.non.zero.entries,rate=0.5) *
                        sample(c(-1,1),num.non.zero.entries,replace=TRUE)
mu[1:14]


###################################################
### code chunk number 6: treethresh.Rnw:178-179
###################################################
x <- mu + rnorm(length(mu))


###################################################
### code chunk number 7: treethresh.Rnw:184-186
###################################################
plot(mu,xlab=paste("Index ",expression(i)),ylab=expression(mu[i]),col=non.zero.entry+1)
title("True signal")


###################################################
### code chunk number 8: treethresh.Rnw:190-192
###################################################
plot(x,xlab=paste("Index ",expression(i)),ylab=expression(x[i]),col=non.zero.entry+1)
title("Observed signal")


###################################################
### code chunk number 9: treethresh.Rnw:200-207
###################################################
values <- c(0,0.01,0.05,0.1,0.2,0.3,0.5,1)
g <- function(x,a=0.5) 0.5*a*exp(a^2/2) * (exp(-a*x) * pnorm(x-a) + exp(a*x)* (1-pnorm(x+a)))
G <- function(x, w=0, a=0.5) (1-w)*pnorm(x)+w*(0.5+integrate(g, 0, x,a=a,rel.tol = .Machine$double.eps^0.25)$value)-0.75

results <- numeric(length(values))
for (i in 1:length(results)) results[i] <- 1/uniroot(G, c(0.5,2), w=values[i], tol = .Machine$double.eps^0.5, maxiter = 2e3)$root



###################################################
### code chunk number 10: treethresh.Rnw:220-222
###################################################
sdev <- mad(x, constant=1.3)
sdev


###################################################
### code chunk number 11: treethresh.Rnw:225-226
###################################################
x <- x /sdev


###################################################
### code chunk number 12: treethresh.Rnw:230-231
###################################################
x.tt <- treethresh(x)


###################################################
### code chunk number 13: treethresh.Rnw:250-251
###################################################
x.tt$splits


###################################################
### code chunk number 14: treethresh.Rnw:257-259
###################################################
# Set the set so that our two calls to prune give the same results
set.seed(111)


###################################################
### code chunk number 15: treethresh.Rnw:262-264
###################################################
x.ttp <- prune(x.tt)
x.ttp$splits


###################################################
### code chunk number 16: treethresh.Rnw:271-274
###################################################
# Set the set so that our two calls to prune give the same results
set.seed(111)
x.ttp <- prune(x.tt)


###################################################
### code chunk number 17: treethresh.Rnw:283-285
###################################################
plot(get.w(x.tt),type="s",xlab=paste("Index ",expression(i)),ylab=expression(hat(w[i]))); abline(v=x.tt$splits[,"pos"],lty=3)
title("Estimated weights (before pruning)")


###################################################
### code chunk number 18: treethresh.Rnw:289-292
###################################################
plot(get.w(x.ttp),type="s",xlab=paste("Index ",expression(i)),ylab=expression(hat(w[i]))); abline(v=x.ttp$splits[,"pos"],lty=3)
title("Estimated weights (after pruning)")
#$


###################################################
### code chunk number 19: treethresh.Rnw:300-301
###################################################
mu.hat <- thresh(x.ttp)


###################################################
### code chunk number 20: treethresh.Rnw:307-309
###################################################
plot(get.t(x.ttp),type="s",xlab=paste("Index ",expression(i)),ylab=expression(t[hat(w[i])])); abline(v=x.ttp$splits[,"pos"],lty=3)
title("Estimated thresholds (after pruning)")


###################################################
### code chunk number 21: treethresh.Rnw:316-317
###################################################
mu.hat <- mu.hat * sdev


###################################################
### code chunk number 22: treethresh.Rnw:324-326
###################################################
plot(mu,xlab=paste("Index ",expression(i)),ylab=expression(hat(mu)[i]),col=non.zero.entry+1)
title("Reconstructed signal")


###################################################
### code chunk number 23: treethresh.Rnw:338-339
###################################################
data(tiles)


###################################################
### code chunk number 24: treethresh.Rnw:344-345
###################################################
set.seed(111)


###################################################
### code chunk number 25: treethresh.Rnw:347-348
###################################################
tiles.noisy <- tiles + 0.8 * rnorm(length(tiles))


###################################################
### code chunk number 26: treethresh.Rnw:355-359
###################################################
zlim <- range(tiles.noisy)
zlim <- c(1.1 * zlim[1] - 0.1 * zlim[2], 1.1 * zlim[2] - 0.1 * zlim[1])
par(mai=rep(0,4))
image(tiles, col=grey(seq(0,1,length.out=256)), zlim=zlim)


###################################################
### code chunk number 27: treethresh.Rnw:363-365
###################################################
par(mai=rep(0,4))
image(tiles.noisy, col=grey(seq(0,1,length.out=256)), zlim=zlim)


###################################################
### code chunk number 28: treethresh.Rnw:373-374
###################################################
tiles.noisy.imwd <- imwd(tiles.noisy)


###################################################
### code chunk number 29: treethresh.Rnw:381-383
###################################################
# Set the seed again to make sure we get the same results as in the next section
set.seed(111)


###################################################
### code chunk number 30: treethresh.Rnw:386-387
###################################################
tiles.noisy.imwd.threshed <- wavelet.treethresh(tiles.noisy.imwd)


###################################################
### code chunk number 31: treethresh.Rnw:394-395
###################################################
tiles.denoised <- imwr(tiles.noisy.imwd.threshed)


###################################################
### code chunk number 32: treethresh.Rnw:401-403
###################################################
par(mai=rep(0,4))
image(tiles.denoised, col=grey(seq(0,1,length.out=256)), zlim=zlim)


###################################################
### code chunk number 33: treethresh.Rnw:407-418
###################################################
sdev <- estimate.sdev(tiles.noisy.imwd)
coefs <- extract.coefficients(tiles.noisy.imwd)
for (nm in names(coefs)) {
  coefs[[nm]] <- coefs[[nm]] / sdev
  coefs[[nm]] <- ebayesthresh(coefs[[nm]], sdev=1)
  coefs[[nm]] <- coefs[[nm]] * sdev
}
tiles.noisy.imwd.threshed.ebs <- insert.coefficients(tiles.noisy.imwd, coefs)
tiles.denoised.ebs <- imwr(tiles.noisy.imwd.threshed.ebs)
par(mai=rep(0,4))
image(tiles.denoised.ebs, col=grey(seq(0,1,length.out=256)), zlim=zlim)


###################################################
### code chunk number 34: treethresh.Rnw:431-432
###################################################
sdev <- estimate.sdev(tiles.noisy.imwd)


###################################################
### code chunk number 35: treethresh.Rnw:437-438
###################################################
tiles.noisy.coefs <- extract.coefficients(tiles.noisy.imwd)


###################################################
### code chunk number 36: treethresh.Rnw:443-445
###################################################
for (nm in names(tiles.noisy.coefs)) 
  tiles.noisy.coefs[[nm]] <- tiles.noisy.coefs[[nm]] / sdev


###################################################
### code chunk number 37: treethresh.Rnw:450-451
###################################################
tiles.noisy.wtt <- wtthresh(tiles.noisy.coefs)


###################################################
### code chunk number 38: treethresh.Rnw:456-458
###################################################
# Set the seed so we can call prune twice and get the same results
set.seed(111)


###################################################
### code chunk number 39: treethresh.Rnw:460-461
###################################################
tiles.noisy.wttp <- prune(tiles.noisy.wtt)


###################################################
### code chunk number 40: treethresh.Rnw:468-470
###################################################
set.seed(111)
tiles.noisy.wttp <- prune(tiles.noisy.wtt)


###################################################
### code chunk number 41: treethresh.Rnw:477-483
###################################################
tiles.noisy.coefs.threshed <- thresh(tiles.noisy.wttp)
for (nm in names(tiles.noisy.coefs)) 
  tiles.noisy.coefs.threshed[[nm]] <-tiles.noisy.coefs.threshed[[nm]] * sdev
tiles.noisy.imwd.threshed <- insert.coefficients(tiles.noisy.imwd, 
                                                   tiles.noisy.coefs.threshed)
tiles.noisy.threshed <- imwr(tiles.noisy.imwd.threshed)


###################################################
### code chunk number 42: treethresh.Rnw:487-507
###################################################
clipped.line <- function(dim, at, clip, factor=1, ...) {
 coords <- clip
 coords[,dim] <- at
 lines(coords*factor+0.5, ...) 
}

draw.cuts <- function(mat, id=1, clip, factor=2, ...) {
  line <- which(mat[,"id"]==id)
  if (is.na(mat[line,"dim"])) 
    return()
  clipped.line(mat[line, "dim"], mat[line, "pos"], clip, factor=factor, ...)
  left.clip <- clip
  left.clip[2,mat[line, "dim"]] <- mat[line,"pos"]
  draw.cuts(mat, mat[line, "left.child.id"], left.clip, factor, ...)
  right.clip <- clip
  right.clip[1,mat[line, "dim"]] <- mat[line,"pos"]
  draw.cuts(mat, mat[line, "right.child.id"], right.clip, factor, ...)  
}

clip <- cbind(c(0,128),c(0,128))


###################################################
### code chunk number 43: treethresh.Rnw:512-519
###################################################
a <- thresh(tiles.noisy.wtt)
for (nm in names(a)) a[[nm]] <- a[[nm]] * sdev
b <- imwr(insert.coefficients(tiles.noisy.imwd, a))
par(mai=rep(0,4))
image(1:128,1:128,b, col=grey(seq(0,1,length.out=256)), zlim=zlim)
invisible(draw.cuts(tiles.noisy.wtt$splits, 1, clip=clip, factor=2, lwd=2))
#$


###################################################
### code chunk number 44: treethresh.Rnw:523-527
###################################################
par(mai=rep(0,4))
image(1:128,1:128,tiles.noisy.threshed, col=grey(seq(0,1,length.out=256)), zlim=zlim)
invisible(draw.cuts(tiles.noisy.wttp$splits, 1, clip=clip, factor=2, lwd=2))
#$


