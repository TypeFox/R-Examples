### R code from vignette source 'rrcov.Rnw'

###################################################
### code chunk number 1: rrcov.Rnw:471-473
###################################################
## set the prompt to "R> " and the continuation to "+ "
options(prompt="R> ", continue="+ ")


###################################################
### code chunk number 2: intro
###################################################
##
## Load the 'rrcov' package and the first two data sets to be
## used throughout the examples
##
library("rrcov")
data("delivery")
delivery.x <- delivery[,1:2]    # take only the X part
data("hbk")
hbk.x <- hbk[,1:3]              # take only the X part


###################################################
### code chunk number 3: intro-mcd
###################################################
##
## Compute MCD estimates for the delivery data set
## - show() and summary() examples
##
mcd <- CovMcd(delivery.x)
mcd
summary(mcd)


###################################################
### code chunk number 4: intro-plot
###################################################
##
## Example plot of the robust against classical
##  distances for the delivery data set
##
plot(mcd, which="dd")


###################################################
### code chunk number 5: intro-S
###################################################
##
## Compute the S-estimates for the delivery data set
##  and provide the standard and extended output
##
est <- CovSest(delivery.x, method="bisquare")
est
summary(est)


###################################################
### code chunk number 6: intro-CovRobust
###################################################
##
## Automatically select the appropriate estimator according
##  to the problem size - in this example the Stahel-Donoho estimates
##  will be selected.
##
est <- CovRobust(delivery.x)
est


###################################################
### code chunk number 7: rrcov.Rnw:1019-1020
###################################################
set.seed(1234)


###################################################
### code chunk number 8: CovControl
###################################################
##
## Controlling the estimation options with a control object
##
control <- CovControlSest(method="biweight")
PcaCov(hbk.x, cov.control=control)



###################################################
### code chunk number 9: CovControl-loop
###################################################
##
## Controlling the estimation options: example
##  of looping through several estimators
##
cc <- list(CovControlMcd(), CovControlMest(), CovControlOgk(), CovControlSest(), CovControlSest(method="rocke"))
clist <- sapply(cc, restimate, x=delivery.x)

sapply(clist, data.class)
sapply(clist, getMeth)


###################################################
### code chunk number 10: CovRobust
###################################################
##
## Automatically select the appropriate estimator according
##  to the problem size.
##
getMeth(CovRobust(matrix(rnorm(40), ncol=2)))       # 20x2      - SDE
getMeth(CovRobust(matrix(rnorm(16000), ncol=8)))    # 2000x8    - bisquare S
getMeth(CovRobust(matrix(rnorm(20000), ncol=10)))   # 2000x10   - Rocke S
getMeth(CovRobust(matrix(rnorm(200000), ncol=2)))   # 100000x2  - OGK


###################################################
### code chunk number 11: CovControl-S
###################################################
##
## Rocke-type S-estimates
##
    getMeth(CovRobust(matrix(rnorm(40), ncol=2), control="rocke"))


###################################################
### code chunk number 12: CovControl-CovRobust
###################################################
##
## Specify some estimation parameters through a control object.
##  The last command line illustrates the accessor method
##  for getting the correlation matrix of the estimate
##  as well as a nice formatting method for covariance
##  matrices.
##
    data("toxicity")
    ctrl <- CovControlOgk(smrob = "s_mad", svrob = "qc")
    est <- CovRobust(toxicity, ctrl)
    round(getCenter(est),2)
    as.dist(round(getCorr(est), 2))


###################################################
### code chunk number 13: ex-cov-plot-lattice
###################################################
##
## Distance plot and Chi-square Q-Q plot of the robust and classical distances
##  Lattice plots (currently not available in rrcov)
##
    library("lattice")
    library("grid")
    library("rrcov")

    data("delivery")
    X <- delivery[,1:2]

    myPanel.dd <- function(x, y, subscripts, cutoff, ...) {
        panel.xyplot(x, y, ...)
        panel.abline(h=cutoff,lty="dashed")

        n <- length(y)
        id.n <- length(which(y>cutoff))
        if(id.n > 0){
            ind <- sort(y, index.return=TRUE)$ix
            ind <- ind[(n-id.n+1):n]

            xrange<-range(y)
            adj <- (xrange[2]-xrange[1])/20
            ltext(x[ind] + adj, y[ind] + adj, ind, cex=1.0)
        }
    }

    myPanel.qchi <- function(x, y, subscripts, cutoff, ...) {
        y <- sort(y, index.return=TRUE)
        iy <- y$ix
        y <- y$x
        panel.xyplot(x, y, ...)
        panel.abline(0,1,lty="dashed")

        n <- length(y)
        id.n <- length(which(y>cutoff))
        if(id.n > 0){
            ind <- (n-id.n+1):n

            xrange<-range(y)
            adj <- (xrange[2]-xrange[1])/50
            ltext(x[ind] + adj, y[ind] + adj, iy[ind], cex=1.0)
        }
    }


    n<-nrow(X)
    p <- ncol(X)
    cutoff <- sqrt(qchisq(0.975, p))

    mm <- covMcd(X)
    dd1 <- sqrt(mm$mah)                             # robust distances

    vv<-cov.wt(X)
    dd2 <- sqrt(mahalanobis(X,vv$center,vv$cov))    # classical distances

    dd<-c(dd1,dd2)                                  # both robust and classical distances

    gr  <- as.factor(c(rep(1,n), rep(2,n)))
    levels(gr)[1] <- "Robust"
    levels(gr)[2] <- "Classical"

    qq  <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

    ind.dd <- c(1:n, 1:n)
    ind.qchi <- c(qq, qq)

    dplot <- xyplot(dd~ind.dd|gr,
                cutoff=cutoff,
                panel = myPanel.dd,
                xlab="Index",
                ylab="Mahalanobis distance",
                main="Distance Plot",
                col = "darkred",
                scales=list(cex=1.0))

    qplot <- xyplot(dd~ind.qchi|gr,
                cutoff=cutoff,
                panel = myPanel.qchi,
                xlab="Quantiles of the chi-squared distribution",
                ylab="Mahalanobis distance",
                main="Chi-Square QQ-Plot",
                col = "darkred",
                scales=list(cex=1.0))

    plot(dplot, split = c(1, 1, 2, 1), more = TRUE)
    plot(qplot, split = c(2, 1, 2, 1), more = FALSE)


###################################################
### code chunk number 14: ex-cov-plot-ellipse
###################################################
##
## a) scatter plot of the data with robust and classical confidence ellipses.
## b) screeplot presenting the robust and classical eigenvalues
##
data("milk")
usr<-par(mfrow=c(1,2))
plot(CovMcd(delivery[,1:2]), which="tolEllipsePlot", classic=TRUE)
plot(CovMcd(milk), which="screeplot", classic=TRUE)
par(usr)


###################################################
### code chunk number 15: pca-ex-hbk
###################################################
##
## Principal Component Analysis example:
##  Plot of the first two principal components of the
##  Hawkins, Bradu and Kass data set: classical and robust
##
par(mfrow=c(1,2))
pca <- PcaClassic(hbk.x, k=2)
cpca <- list(center=c(0,0), cov=diag(pca@eigenvalues), n.obs=pca@n.obs)
rrcov:::.myellipse(pca@scores, xcov=cpca, main="Classical",  xlab="PC1", ylab="PC2", ylim=c(-30,30), id.n=0)
abline(v=0)
abline(h=0)
text(-29,6,labels="1-10", cex=0.8)
text(-37,-13,labels="14", cex=0.8)
text(-31,-5,labels="11-13", cex=0.8)
hub <- PcaHubert(hbk.x, k=2, mcd=TRUE)
chub <- list(center=c(0,0), cov=diag(hub@eigenvalues), n.obs=hub@n.obs)
rrcov:::.myellipse(hub@scores, xcov=chub, main="Robust (MCD)", xlab="PC1", ylab="PC2", xlim=c(-10,45), ylim=c(-25,15), id.n=0)
abline(v=0)
abline(h=0)
text(30,-9,labels="1-10", cex=0.8)
text(36,-11,labels="11", cex=0.8)
text(42,-4,labels="12", cex=0.8)
text(41,-11,labels="13", cex=0.8)
text(44,-15,labels="14", cex=0.8)


###################################################
### code chunk number 16: pca-classic
###################################################
##
## Classical PCA
##
pca <- PcaClassic(~., data=hbk.x)
pca
summary(pca)
plot(pca)
getLoadings(pca)


###################################################
### code chunk number 17: pca-robust
###################################################
##
## Robust PCA
##
rpca <- PcaGrid(~., data=hbk.x)
rpca
summary(rpca)


###################################################
### code chunk number 18: ex-pca-plot-screeplot
###################################################
##
## Screeplot for classical and robust PCA of the milk data set.
##
usr <- par(mfrow=c(1,2))
screeplot(PcaClassic(milk), type="lines",
main="Screeplot: classical PCA", sub="milk data")
screeplot(PcaHubert(milk), type="lines", main="Screeplot: robust PCA",
sub="milk data")
par(usr)


###################################################
### code chunk number 19: ex-pca-plot-biplot
###################################################
##
## Robust biplot for the UN86 data
##
data("un86")
set.seed(9)

usr<-par(mfrow=c(1,2))
biplot(PcaCov(un86, corr=TRUE, cov.control=NULL),
main="Classical biplot", col=c("gray55", "red"))
biplot(PcaCov(un86, corr=TRUE), main="Robust biplot",
col=c("gray55", "red"))
par(usr)


###################################################
### code chunk number 20: ex-pca-plot-diagplot
###################################################
##
## An example of the classical and robust diagnostic
##  plot for the hbk data set
##
usr<-par(mfrow=c(1,2))
plot(PcaClassic(hbk.x, k=2), sub="data set: hbk, k=2")
plot(PcaHubert(hbk.x, k=2), sub="data set: hbk, k=2")
par(usr)


###################################################
### code chunk number 21: ex-pca-plot-distplot
###################################################
##
## If k=p the orthogonal distances are not meaningful and
##  simple distance plot of the score distances will be shown
##
usr<-par(mfrow=c(1,2))
plot(PcaClassic(hbk.x, k=3), sub="data set: hbk.x, k=3")
plot(PcaHubert(hbk.x, k=3), sub="data set: hbk.x, k=3")
par(usr)


###################################################
### code chunk number 22: lda-prelim
###################################################
data("diabetes", package="mclust")


###################################################
### code chunk number 23: lda-cloud
###################################################
library("lattice")    # load the graphical library

## set different plot symbols - important for black-and-white print
sup.sym <- trellis.par.get("superpose.symbol")
sup.sym$pch <- c(1,2,3,4,5)
trellis.par.set("superpose.symbol", sup.sym)
cloud.plt <- cloud(insulin ~ glucose + sspg, groups = class, data = diabetes, auto.key=TRUE)


###################################################
### code chunk number 24: lda-cloud-fig
###################################################
print(cloud.plt)


###################################################
### code chunk number 25: lda-classic
###################################################
lda <- LdaClassic(class~., data=diabetes)
lda


###################################################
### code chunk number 26: lda-classic-predict
###################################################
predict(lda)


###################################################
### code chunk number 27: lda-robust
###################################################
rlda <- Linda(class~., data=diabetes)
rlda

rlda.predict <- predict(rlda)
cat("\nApparent error rate: ", round(rrcov:::.AER(rlda.predict@ct),4))


###################################################
### code chunk number 28: qda-robust
###################################################
rqda <- QdaCov(class~., data=diabetes)
rqda
rqda.predict <- predict(rqda)
cat("\nApparent error rate: ", round(rrcov:::.AER(rqda.predict@ct),4))


