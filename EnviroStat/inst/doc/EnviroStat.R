### R code from vignette source 'EnviroStat.Rnw'

###################################################
### code chunk number 1: EnviroStat.Rnw:68-70
###################################################
library(EnviroStat)
data(ozone.NY)


###################################################
### code chunk number 2: EnviroStat.Rnw:75-78
###################################################
month <- ozone.NY[,1]
weekday <- ozone.NY[,2]
sqO3 <- ozone.NY[,3:38]


###################################################
### code chunk number 3: EnviroStat.Rnw:84-85
###################################################
data(location.NY)


###################################################
### code chunk number 4: EnviroStat.Rnw:93-94
###################################################
sitenumber <- apply(is.na(sqO3), 2, sum)


###################################################
### code chunk number 5: EnviroStat.Rnw:100-105
###################################################
norder <- c(2, 6, 5, 1, 3, 4, 7, 8, 9)
tt <- NULL
for (i in 1:9) tt <- c(tt, c(1:4) + 4 * (norder[i]-1))
ndata <- sqO3[,tt]
nloc <- location.NY[norder, ]


###################################################
### code chunk number 6: EnviroStat.Rnw:110-113
###################################################
plot(nloc[,2], nloc[,1], type = "n",
     xlab = "Long", ylab = "Lat")
text(nloc[,2], nloc[,1], c(1:9))


###################################################
### code chunk number 7: EnviroStat.Rnw:143-145
###################################################
month[1:5]
weekday[1:5]


###################################################
### code chunk number 8: EnviroStat.Rnw:151-152
###################################################
ZZ <- model.matrix(~as.factor(month) + as.factor(weekday))


###################################################
### code chunk number 9: EnviroStat.Rnw:226-228
###################################################
em.fit <- staircase.EM(ndata, p = 4, covariate = ZZ,
                       maxit = 200, tol = .000001)


###################################################
### code chunk number 10: EnviroStat.Rnw:232-233
###################################################
em.fit$block


###################################################
### code chunk number 11: EnviroStat.Rnw:238-239
###################################################
em.fit$Omega


###################################################
### code chunk number 12: EnviroStat.Rnw:245-246 (eval = FALSE)
###################################################
## em.fit$Lambda


###################################################
### code chunk number 13: EnviroStat.Rnw:260-261 (eval = FALSE)
###################################################
## em.fit$Beta0[,1:4]


###################################################
### code chunk number 14: EnviroStat.Rnw:267-273
###################################################
cov.est <- em.fit$Psi[[1]]
dim1 <- nrow(cov.est)
dim2 <- nrow(em.fit$Omega)
corr.est <- cov.est / 
    sqrt( matrix(diag(cov.est), dim1, dim1) *
          t(matrix(diag(cov.est), dim1, dim1)) )


###################################################
### code chunk number 15: EnviroStat.Rnw:278-279 (eval = FALSE)
###################################################
## round(corr.est, 2)


###################################################
### code chunk number 16: EnviroStat.Rnw:286-290
###################################################
omega <- em.fit$Omega / 
    sqrt( matrix(diag(em.fit$Omega), dim2, dim2) *
          t(matrix(diag(em.fit$Omega), dim2, dim2)) )
round(omega, 2)


###################################################
### code chunk number 17: EnviroStat.Rnw:307-308
###################################################
coords <- Flamb2(abs(nloc))


###################################################
### code chunk number 18: EnviroStat.Rnw:314-315
###################################################
dist <- Fdist(coords$xy)


###################################################
### code chunk number 19: EnviroStat.Rnw:320-328
###################################################
par(mfrow = c(1, 1))
plot(-.2, 0, xlim = c(0, 250), ylim = c(-.2, 1),
     xlab = "Dist", ylab = "Spatial correlation", type = "n",
     main = c("Intersite correlation vs intersite distance",
            "(via the Lambert projection)"))
for (i in 1:8) 
    for (j in (i+1):9)
        points(dist[i, j], corr.est[i, j])                         


###################################################
### code chunk number 20: dispersion-dist
###################################################
disp <- 2 - 2 * corr.est
plot(-.2, 0, xlim = c(0, 250), ylim = c(0, 2),
     xlab = "Dist", ylab = "Dispersion", type = "n",
     main = c("Intersite dispersion vs intersite distance",
            "(via the Lambert projection)"))
for (i in 1:8)
    for (j in (i+1):9)
        points(dist[i, j], disp[i, j])


###################################################
### code chunk number 21: EnviroStat.Rnw:365-373
###################################################
h.lt <- dist[row(dist) < col(dist)]
disp.lt <- disp[row(disp) < col(disp)]	
variogfit <- Fvariogfit3(disp.lt, h.lt, a0 = 1.5, t0 = .1)
x <- seq(min(h.lt), max(h.lt), 1)
a0 <- variogfit$a[1]
t0 <- variogfit$t0
disp <- 2 - 2 * corr.est
plot(-.2, 0, xlim = c(0, 250), ylim = c(0, 2),
     xlab = "Dist", ylab = "Dispersion", type = "n",
     main = c("Intersite dispersion vs intersite distance",
            "(via the Lambert projection)"))
for (i in 1:8)
    for (j in (i+1):9)
        points(dist[i, j], disp[i, j])
lines(x, a0 + (2 - a0) * (1 - exp(-(t0 * x))))


###################################################
### code chunk number 22: EnviroStat.Rnw:442-443
###################################################
coords.lamb  <- coords$xy / 10 


###################################################
### code chunk number 23: EnviroStat.Rnw:447-449
###################################################
sg.est <- Falternate3(disp, coords.lamb, max.iter = 100,
                      alter.lim = 100, model = 1)


###################################################
### code chunk number 24: EnviroStat.Rnw:467-475 (eval = FALSE)
###################################################
## apply(coords.lamb, 2, range)
## coords.grid <- Fmgrid(range(coords.lamb[,1]),
##                       range(coords.lamb[,2]))
## par(mfrow = c(1, 2))
## temp <- setplot(coords.lamb, axis = TRUE)  
## deform  <- Ftransdraw(disp = disp, Gcrds = coords.lamb,
##                       MDScrds = sg.est$ncoords,
##                       gridstr = coords.grid)


###################################################
### code chunk number 25: EnviroStat.Rnw:493-494
###################################################
Tspline <- sinterp(coords.lamb, sg.est$ncoords, lam = 50 )


###################################################
### code chunk number 26: EnviroStat.Rnw:508-514 (eval = FALSE)
###################################################
## par(mfrow = c(1, 1))
## Tgrid  <- bgrid(start = c(0, 0), xmat = coords.lamb,
##                 coef = Tspline$sol)
## tempplot <- setplot(coords.lamb, axis = TRUE)
## text(coords.lamb, labels = 1:nrow(coords.lamb))
## draw(Tgrid, fs = TRUE)


###################################################
### code chunk number 27: EnviroStat.Rnw:527-530
###################################################
lat10 <- seq(min(nloc[,1]), max(nloc[,1]), length = 10)
long10 <- seq(max(abs(nloc[,2])), min(abs(nloc[,2])), length = 10)
llgrid <- cbind(rep(lat10, 10), c(outer(rep(1, 10), long10)))


###################################################
### code chunk number 28: EnviroStat.Rnw:542-545
###################################################
z <- coords
newcrds.lamb <- Flamb2(llgrid, latrf1 = z$latrf1, latrf2 = z$latrf2,
                       latref = z$latref, lngref = z$lngref)$xy / 10


###################################################
### code chunk number 29: EnviroStat.Rnw:551-552
###################################################
allcrds <- rbind(newcrds.lamb, coords.lamb)


###################################################
### code chunk number 30: EnviroStat.Rnw:557-560
###################################################
corr.est <- corrfit(allcrds, Tspline = Tspline,
                    sg.fit  = sg.est, model = 1)
round(corr.est$cor[1:5, 1:5], 2)


###################################################
### code chunk number 31: EnviroStat.Rnw:574-575
###################################################
diag(cov.est)


###################################################
### code chunk number 32: EnviroStat.Rnw:580-583
###################################################
Tspline.var <- sinterp(allcrds[101:109,],
                       matrix(diag(cov.est), ncol = 1),
                       lam = 50)


###################################################
### code chunk number 33: EnviroStat.Rnw:589-591
###################################################
varfit <- seval(allcrds, Tspline.var)$y
temp <- matrix(varfit, length(varfit), length(varfit))


###################################################
### code chunk number 34: EnviroStat.Rnw:597-598
###################################################
covfit <- corr.est$cor * sqrt(temp * t(temp))


###################################################
### code chunk number 35: EnviroStat.Rnw:611-615
###################################################
u <- 100 # number of new locations
p <- 4  # dimension of the multivariate response
hyper.est <- staircase.hyper.est(emfit = em.fit,
                                 covfit = covfit, u = u, p = p)


###################################################
### code chunk number 36: EnviroStat.Rnw:623-627
###################################################
x <- hyper.est
tpt <- 183
Z <- x$covariate[tpt,]
y <- x$data[tpt,]


###################################################
### code chunk number 37: EnviroStat.Rnw:632-633
###################################################
b0 <- matrix(rep(c(x$Beta0[,1:4]), 100), nrow = 12)


###################################################
### code chunk number 38: EnviroStat.Rnw:637-639
###################################################
mu.u <- Z %*% b0 + (as.matrix(y)- Z %*% x$Beta0) %*% 
                    kronecker(x$Xi0.0, diag(4))


###################################################
### code chunk number 39: EnviroStat.Rnw:644-646
###################################################
color <- colors()[c(12, 26, 32, 37, 53, 60, 70, 80, 84, 88, 94, 101,
                    116, 142, 366, 371, 376, 386, 392, 398, 400:657)]


###################################################
### code chunk number 40: EnviroStat.Rnw:650-655
###################################################
par(mfrow = c(1, 1))
plot(c(7, 13), range(y), type = "n", 
     xlab = "Hours", ylab = "Levels (log)")
for (i in 1:4)
    points(rep(i+7, 9), y[i + 4 * 0:8], col = color)


###################################################
### code chunk number 41: EnviroStat.Rnw:660-661
###################################################
pdf('EnviroStat-41-foo.pdf')


###################################################
### code chunk number 42: EnviroStat.Rnw:664-673
###################################################
par(mfrow = c(2, 2))
for (i in 1:4) {
    tt <- i + 4 * 0:99
    mu <- mu.u[tt]
    hr <- matrix(mu, byrow = TRUE, ncol = 10)
    print(range(hr))
    contour(-long10, lat10, hr, xlab = "Long", ylab = "Lat", 
            main = paste("Mean: Day 183 - Hour ", 7+i))
}


###################################################
### code chunk number 43: EnviroStat.Rnw:676-677
###################################################
dev.off()


###################################################
### code chunk number 44: EnviroStat.Rnw:699-700
###################################################
simu <- pred.dist.simul(hyper.est, tpt = 183, N = 1000)


###################################################
### code chunk number 45: EnviroStat.Rnw:705-706
###################################################
x <- apply(simu, 2, mean)[1:400]


###################################################
### code chunk number 46: EnviroStat.Rnw:710-711
###################################################
pdf('EnviroStat-44-foo.pdf')


###################################################
### code chunk number 47: EnviroStat.Rnw:714-723
###################################################
par(mfrow = c(2, 2))
for (i in 1:4) {
    tt <- i + 4 * 0:99
    x1 <- x[tt]
    hr <- matrix(x1, byrow = TRUE, ncol = 10)
    print(range(x1 ))
    contour(-long10, lat10, hr, xlab = "Long", ylab = "Lat", 
            main = paste("Mean: Day 183 - Hour ", 7+i))
}


###################################################
### code chunk number 48: EnviroStat.Rnw:726-727
###################################################
dev.off()


###################################################
### code chunk number 49: EnviroStat.Rnw:734-735
###################################################
x <- simu[,1:400]


###################################################
### code chunk number 50: EnviroStat.Rnw:739-740
###################################################
pdf('EnviroStat-46-foo.pdf')


###################################################
### code chunk number 51: EnviroStat.Rnw:743-753
###################################################
par(mfrow = c(2, 2))
for (i in 1:4) {
    tt <- i + 4 * 0:99
    x1 <- x[,tt]
    x2 <- diag(var(x1))
    vv <- matrix(x2, byrow = TRUE, ncol = 10)
    contour(-long10, lat10, vv, xlab = "Long", ylab = "Lat", 
            main = paste("Var: Day 183 - Hour ", 7+i))
points(nloc[,2], nloc[,1])
}


###################################################
### code chunk number 52: EnviroStat.Rnw:756-757
###################################################
dev.off()


###################################################
### code chunk number 53: EnviroStat.Rnw:784-787
###################################################
nsel <- 3
yy <- ldet.eval((hyper.est$Lambda.0 + t(hyper.est$Lambda.0)) / 2,
                nsel, all = FALSE)


###################################################
### code chunk number 54: EnviroStat.Rnw:793-796 (eval = FALSE)
###################################################
## yy1 <- ldet.eval(((hyper.est$Lambda.0 + 
##                    t(hyper.est$Lambda.0))/2)[1:10, 1:10],
##                  nsel, all = TRUE)


