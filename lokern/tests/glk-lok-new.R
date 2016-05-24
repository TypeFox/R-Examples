require(lokern)
data(xSim)
(n <- length(xSim))
summary(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x

lk  <- lokerns(tt, xSim, trace=1)
head(x. <- seq(-3.5, 4.5, length=1001))
lk.  <- lokerns(tt, xSim, x.out= x., trace=1)
lkw <- tryCatch(lokerns(tt, xSim, x.out= x., x.inOut="aim"),
                warning = function(.) .)
stopifnot(inherits(lkw, "simpleWarning"),
	  grepl("x.inOut", lkw[["message"]]))

lkn. <- lokerns(tt, xSim, n.out= 1000, x.inOut=FALSE)
lknT <- lokerns(tt, xSim, n.out= 1000, x.inOut=TRUE)
lkn1 <- lokerns(tt, xSim, n.out= 1000, x.inOut="aim")
lkn2 <- lokerns(tt, xSim, n.out= 1000, x.inOut="simple")
lkn3 <- lokerns(tt, xSim, n.out= 1000, x.inOut="interpolate")

iO <- which("call" != names(lkn.))
stopifnot(all.equal(lknT[iO],
		    lkn1[iO], tol=1e-12),
          inherits(try(residuals(lkn.), silent=TRUE), "try-error"),
          length(residuals(lkn1)) == length(tt),
          length(fitted   (lkn2)) == length(tt),
          TRUE)
pdf("lok-xout-ex.pdf")

plot(tt, xSim, main = "data + 4 different lokerns(.., x.inOut = *)")
lines(est~x.out, data=lk. [c("x.out", "est")])
lines(est~x.out, data=lkn1[c("x.out", "est")], col=2)
lines(est~x.out, data=lkn2[c("x.out", "est")], col=3)
lines(est~x.out, data=lkn3[c("x.out", "est")], col=4)
mtext("(you should not see different curves -- rather perfectly overlap")

gkn. <- glkerns(tt, xSim, n.out= 1000, x.inOut=FALSE)
gknT <- glkerns(tt, xSim, n.out= 1000, x.inOut=TRUE)
gkn1 <- glkerns(tt, xSim, n.out= 1000, x.inOut="aim")
gkn2 <- glkerns(tt, xSim, n.out= 1000, x.inOut="simple")
gkn3 <- glkerns(tt, xSim, n.out= 1000, x.inOut="interpolate")

plot(tt, xSim, main = "data + 4 different glkerns(.., x.inOut = *)")
lines(est~x.out, data=gkn.[c("x.out", "est")])
lines(est~x.out, data=gkn1[c("x.out", "est")], col=2)
lines(est~x.out, data=gkn2[c("x.out", "est")], col=3)
lines(est~x.out, data=gkn3[c("x.out", "est")], col=4)
mtext("(you should not see different curves -- rather perfectly overlap")

stopifnot(all(sapply(ls(patt="kn"), function(.)get(.)$sig) -
	      0.145106481 < 5e-10))

cat('Time elapsed: ', proc.time(),'\n') # "stats"
