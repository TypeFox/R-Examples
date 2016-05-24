### R code from vignette source 'allometry.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(width=65)
options(repos="http://cran.r-project.org")
#if(!require(systemfit, quietly=TRUE)) 
#  install.packages("systemfit", dependencies=TRUE)
#if(!require(nlreg, quietly=TRUE)) 
#  install.packages("nlreg", dependencies=TRUE)
require(MASS, quietly=TRUE)
require(lattice, quietly=TRUE)
lattice.options(default.theme = canonical.theme(color = FALSE))
require(nlme, quietly=TRUE)
require(mgcv, quietly=TRUE)


###################################################
### code chunk number 2: fig-norm-dd
###################################################
set.seed(1)
diameters <- rnorm(30, mean = 50, sd = 6)
hist(diameters, main = "", xlab = "Diameter Class (cm)")


###################################################
### code chunk number 3: norm-dd
###################################################
par(las=1, mar=c(5,5,1,2))
x.bar <- mean(diameters)
s.x <- sd(diameters)
set.seed(1)
diameters <- rnorm(30, mean = 50, sd = 6)
hist(diameters, main = "", xlab = "Diameter Class (cm)")
curve(dnorm(x, 50, 6) * 30 * 6, add=TRUE, col="darkgrey")
curve(dnorm(x, x.bar, s.x) * 30 * s.x, add=TRUE, lty=3)
legend("topleft", lty = c(1,3), bty = "n", cex=0.8, 
col=c("darkgrey","black"), legend = c("True pdf","Fitted pdf"))


###################################################
### code chunk number 4: allometry.rnw:188-191
###################################################
min(diameters)
quantile(diameters, p = (5 + (0:9) * 10) / 100 )
max(diameters)


###################################################
### code chunk number 5: allometry.rnw:253-255
###################################################
(x.bar.pp <- mean(diameters))
(s.x.pp <- sd(diameters))


###################################################
### code chunk number 6: allometry.rnw:270-273
###################################################
curve(dnorm(x, x.bar.pp, s.x.pp) * length(diameters) * s.x.pp, 
      add = TRUE, 
      col = "darkgrey")


###################################################
### code chunk number 7: allometry.rnw:342-345
###################################################
(x.bar.mm <- mean(diameters))
(s.x.mm <- sqrt(sum((diameters - x.bar.mm)^2) / 
                length(diameters)))


###################################################
### code chunk number 8: allometry.rnw:375-380
###################################################
ll.norm <- function(parameters, data) {
  mu.hat <- parameters[1]
  sigma.hat <- parameters[2]
  return(sum(dnorm(data, mu.hat, sigma.hat, log = TRUE)))
}


###################################################
### code chunk number 9: allometry.rnw:391-396
###################################################
(mle.n <- optim(c(x.bar.mm = x.bar.mm, s.x.mm = s.x.mm), 
             ll.norm, 
             data = diameters, 
	     hessian = TRUE,
             control = list(fnscale = -1)))


###################################################
### code chunk number 10: allometry.rnw:415-416
###################################################
mle.n$par


###################################################
### code chunk number 11: allometry.rnw:423-424
###################################################
solve(-mle.n$hessian)


###################################################
### code chunk number 12: allometry.rnw:429-430
###################################################
sqrt(diag(solve(-mle.n$hessian)))


###################################################
### code chunk number 13: allometry.rnw:456-459
###################################################
library(MASS)
d.n <- fitdistr(diameters, "normal")
d.n


###################################################
### code chunk number 14: allometry.rnw:479-481
###################################################
d.w2 <- fitdistr(diameters, "weibull")
d.w2


###################################################
### code chunk number 15: allometry.rnw:490-492
###################################################
logLik(d.n)
logLik(d.w2)


###################################################
### code chunk number 16: allometry.rnw:513-517
###################################################
dweibull3 <- function(x, gamma, beta, alpha) {
  (gamma/beta)*((x - alpha)/beta)^(gamma - 1) * 
    (exp(-((x - alpha)/beta)^gamma))
}


###################################################
### code chunk number 17: allometry.rnw:523-525
###################################################
ll.w3 <- function(p, data) 
     sum(log(dweibull3(data, p[1], p[2], p[3])))


###################################################
### code chunk number 18: allometry.rnw:567-573
###################################################
mle.w3.nm <- optim(c(gamma = 1, beta = 5, alpha = 10), 
                   ll.w3, 
                   data = diameters,
                   hessian = TRUE,
                   control = list(fnscale = -1))
mle.w3.nm$par


###################################################
### code chunk number 19: allometry.rnw:577-578
###################################################
mle.w3.nm$value


###################################################
### code chunk number 20: allometry.rnw:584-585
###################################################
sqrt(diag(solve(-mle.w3.nm$hessian)))


###################################################
### code chunk number 21: allometry.rnw:596-606
###################################################
mle.w3.bfgs <- optim( c(gamma = 1, beta = 5, alpha = 10), 
                      ll.w3, 
                      method = "BFGS",
                      data = diameters,
                      hessian = TRUE,
                      control = list(fnscale = -1, 
                                     maxit = 1000))
mle.w3.bfgs$value
mle.w3.bfgs$par
(mle.w3.bfgs.se <- sqrt(diag(solve(-mle.w3.bfgs$hessian))))


###################################################
### code chunk number 22: fig-3w-dd
###################################################
hist(diameters, main="", xlab="Diam. Class (cm)", freq=FALSE)
curve(dnorm(x, x.bar, s.x), add=TRUE, col=grey(0.2), lty=1)
curve(dweibull(x, d.w2$estimate[1], d.w2$estimate[2]), 
      add=TRUE, col=grey(0.2), lty = 2)
w3.n.h <- mle.w3.nm$par
curve(dweibull(x - w3.n.h[3], w3.n.h[1], w3.n.h[2]), 
      add=TRUE, col=grey(0.2), lty = 3)
w3.b.h <- mle.w3.bfgs$par
curve(dweibull(x - w3.b.h[3], w3.b.h[1], w3.b.h[2]), 
      add=TRUE, col=grey(0.2), lty = 4)
legend("topleft", lty = 1:4, bty = "n", cex=0.8,
       legend = c("Normal","Weibull (2)",
         "W. (3, NM)","W. (3, BFGS)")) 


###################################################
### code chunk number 23: 3w-dd
###################################################
par(las=1, mar=c(5,5,1,0))
hist(diameters, main="", xlab="Diam. Class (cm)", freq=FALSE)
curve(dnorm(x, x.bar, s.x), add=TRUE, col=grey(0.2), lty=1)
curve(dweibull(x, d.w2$estimate[1], d.w2$estimate[2]), 
      add=TRUE, col=grey(0.2), lty = 2)
w3.n.h <- mle.w3.nm$par
curve(dweibull(x - w3.n.h[3], w3.n.h[1], w3.n.h[2]), 
      add=TRUE, col=grey(0.2), lty = 3)
w3.b.h <- mle.w3.bfgs$par
curve(dweibull(x - w3.b.h[3], w3.b.h[1], w3.b.h[2]), 
      add=TRUE, col=grey(0.2), lty = 4)
legend("topleft", lty = 1:4, bty = "n", cex=0.8,
       legend = c("Normal","Weibull (2)",
         "W. (3, NM)","W. (3, BFGS)")) 


###################################################
### code chunk number 24: allometry.rnw:676-688
###################################################
all.w3 <- list(loc = function(p, data, fix) 
               sum(log(dweibull3(data, p[1], p[2], fix))),
     	       sha = function(p, data, fix) 
               sum(log(dweibull3(data, fix, p[1], p[2]))),
     	       sca = function(p, data, fix) 
               sum(log(dweibull3(data, p[1], fix, p[2]))),
               loc.sha = function(p, data, fix) 
               sum(log(dweibull3(data, fix[1], p, fix[2]))),
     	       sca.sha = function(p, data, fix) 
               sum(log(dweibull3(data, fix[1], fix[2], p))),
     	       sca.loc = function(p, data, fix) 
               sum(log(dweibull3(data, p, fix[1], fix[2]))))


###################################################
### code chunk number 25: allometry.rnw:695-710
###################################################
grain <- 30
k <- 1

p.hat <- mle.w3.bfgs$par
p.se <- mle.w3.bfgs.se

frame <- list(sha = seq(from = p.hat[1] - k * p.se[1], 
                 to = p.hat[1] + k * p.se[1],
                 length = grain ),
              sca = seq(from = p.hat[2] - k * p.se[2], 
                to = p.hat[2] + k * p.se[2],
                length = grain ),
              loc = seq(from = p.hat[3] - k * p.se[3], 
                to = p.hat[3] + k * p.se[3],
                length = grain ))


###################################################
### code chunk number 26: allometry.rnw:726-737
###################################################
profile.fn <- function(x, fn, p.fix) {
  out <- try(optim( p.hat[-p.fix], 
                   fn, 
                   method = "BFGS",
                   data = diameters,
                   fix = x,
                   control = list(fnscale = -1, 
                     maxit = 1000)), silent = TRUE)
  if (class(out) == "try-error") return(NA)
  else return(out$value)
}


###################################################
### code chunk number 27: allometry.rnw:743-755
###################################################
sha.out <- 
  sapply(frame$sha, profile.fn, fn = all.w3$sha, p.fix = 1)
sca.out <- 
  sapply(frame$sca, profile.fn, fn = all.w3$sca, p.fix = 2)
loc.out <- 
  sapply(frame$loc, profile.fn, fn = all.w3$loc, p.fix = 3)
ll.stack <- as.data.frame(rbind(cbind(sha.out, frame$sha),
                                cbind(sca.out, frame$sca),
                                cbind(loc.out, frame$loc)))               
names(ll.stack) <- c("ll","x")
ll.stack$parameter <- rep(c("shape","scale","location"), 
                          rep(grain,3))


###################################################
### code chunk number 28: 3.profiles
###################################################
xyplot(ll ~ x | parameter, 
             type = "l", 
             ylab = "ll(x)",
             scales = list(x="free"),
             layout = c(3,1), 
             data = ll.stack)


###################################################
### code chunk number 29: profile-1
###################################################
print(
xyplot(ll ~ x | parameter, 
             type = "l", 
             ylab = "ll(x)",
             scales = list(x="free"),
             layout = c(3,1), 
             data = ll.stack)
      )


###################################################
### code chunk number 30: allometry.rnw:802-803
###################################################
sha.sca.grid <- expand.grid(sha = frame$sha, sca = frame$sca)


###################################################
### code chunk number 31: fit
###################################################
sha.sca.grid$ll <- 
  with(sha.sca.grid, 
       mapply(function(sha, sca) {
         out <- try(optim(p.hat[3], 
                          all.w3$sca.sha, 
                          method = "BFGS",
                          data = diameters,
                          fix = c(sha, sca),
                          control = list(fnscale = -1, 
                            maxit = 1000)))
         if (class(out) == "try-error") return(NA)
         else return(out$value)
       }, sha, sca))


###################################################
### code chunk number 32: allometry.rnw:827-828
###################################################
sha.sca.grid <- na.omit(sha.sca.grid)


###################################################
### code chunk number 33: fig-profile-2
###################################################
contourplot(ll ~ sha * sca, data = sha.sca.grid,
            at = c(-92, -96,-100, -110),
            xlab = "Shape", ylab = "Scale",
            panel = function(x, y, z, ...) {  
              panel.contourplot(x, y, z, ...)
              panel.points(x, y, col="darkgrey")
              panel.points(p.hat[1], p.hat[2], 
	                   pch = 19, col = "black")
            })


###################################################
### code chunk number 34: profile-2
###################################################
print(
contourplot(ll ~ sha * sca, data = sha.sca.grid,
            at = c(-92, -96,-100, -110),
            xlab = "Shape", ylab = "Scale",
            panel = function(x, y, z, ...) {  
              panel.contourplot(x, y, z, ...)
              panel.points(x, y, col="darkgrey")
              panel.points(p.hat[1], p.hat[2], 
	                   pch = 19, col = "black")
            })
)


###################################################
### code chunk number 35: allometry.rnw:910-919
###################################################
ll.w3.sb2 <- function(p, data) {
  gam <- p[1]
  bet <- p[2]
  alp <- p[3]
  mu.prime.2 <- bet^2 * gamma(2/gam + 1) + 
    2 * bet * gamma(1/gam + 1) * alp + alp^2
  sum(log(dweibull3(data, gam, bet, alp))) - 
    log(mu.prime.2) * length(data)
}


###################################################
### code chunk number 36: allometry.rnw:926-933
###################################################
mle.w3.sb2 <- optim( c(gamma = 1, beta = 5, alpha = 0), 
                     ll.w3.sb2, 
                     method = "BFGS",
                     data = diameters,
                     hessian = TRUE,
                     control = list(fnscale = -1, 
                                    maxit = 1000))


###################################################
### code chunk number 37: allometry.rnw:937-939
###################################################
mle.w3.sb2$par
sqrt(diag(solve(-mle.w3.sb2$hessian)))


###################################################
### code chunk number 38: another-try (eval = FALSE)
###################################################
## set.seed(1)
## diameters <- rnorm(300, mean = 50, sd = 6)
## 
## mle.w3.bfgs <- optim( c(1, 5, 10), 
##                      ll.w3, 
##                      method = "BFGS",
##                      data = diameters,
##                      hessian = TRUE,
##                      control = list(fnscale = -1, 
##                        maxit = 1000))
## mle.w3.bfgs.se <- sqrt(diag(solve(-mle.w3.bfgs$hessian)))
## 
## p.hat <- mle.w3.bfgs$par
## p.se <- mle.w3.bfgs.se
## 
## k <- 2
## 
## frame <- list( sha = seq(from = p.hat[1] - k * p.se[1], 
##                  to = p.hat[1] + k * p.se[1] ,
##                  length = grain ),
##               sca = seq(from = p.hat[2] - k * p.se[2], 
##                 to = p.hat[2] + k * p.se[2] ,
##                 length = grain ),
##               loc = seq(from = p.hat[3] - k * p.se[3], 
##                 to = p.hat[3] + k * p.se[3] ,
##                 length = grain ))
## 
## sha.sca.grid <- expand.grid(sha = frame$sha, sca = frame$sca)
## 
## sha.sca.grid <- expand.grid(sha = frame$sha, sca = frame$sca)
## 
## sha.sca.grid$ll <- 
##   with(sha.sca.grid, 
##        mapply(function(sha, sca) {
##          out <- try(optim( p.hat[3], 
##                           all.w3$sca.sha, 
##                           method = "BFGS",
##                           data = diameters,
##                           fix = c(sha, sca),
##                           control = list(fnscale = -1, 
##                             maxit = 1000)), silent=FALSE)
##          if (class(out) == "try-error") return(NA)
##          else return(out$value)
##        }, sha, sca))
## 
## sha.sca.grid <- na.omit(sha.sca.grid)


###################################################
### code chunk number 39: fig-profile-3 (eval = FALSE)
###################################################
## print(contourplot(ll ~ sha * sca, data = sha.sca.grid,
##                   at = max(sha.sca.grid$ll) - (0:4) * 4,
##                   panel = function(x,y,z,...) {
## 	    	    panel.contourplot(x,y,z,...)
##  		    panel.points(x, y, col="darkgrey")
##                     panel.points(p.hat[1],p.hat[2], pch=19, col="black")
##                   }))
## 


###################################################
### code chunk number 40: profile-3 (eval = FALSE)
###################################################
## print(contourplot(ll ~ sha * sca, data = sha.sca.grid,
##                   at = max(sha.sca.grid$ll) - (0:4) * 4,
##                   panel = function(x,y,z,...) {
## 	    	    panel.contourplot(x,y,z,...)
##  		    panel.points(x, y, col="darkgrey")
##                     panel.points(p.hat[1],p.hat[2], pch=19, col="black")
##                   }))
## 


###################################################
### code chunk number 41: allometry.rnw:1026-1028
###################################################
system("rm -fr package-Ch5")
package.skeleton(name = "package-Ch5")


