### R code from vignette source 'selection.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: selection.Rnw:67-68
###################################################
options( prompt = "R> ", ctinue = "+  " )


###################################################
### code chunk number 2: code:t2generate
###################################################
  set.seed(0)
  library("sampleSelection")
  library("mvtnorm")
  eps <- rmvnorm(500, c(0,0), matrix(c(1,-0.7,-0.7,1), 2, 2))
  xs <- runif(500)
  ys <- xs + eps[,1] > 0
  xo <- runif(500)
yoX <- xo + eps[,2]
yo <- yoX*(ys > 0)


###################################################
### code chunk number 3: code:t2summary
###################################################
  summary( selection(ys~xs, yo ~xo))


###################################################
### code chunk number 4: selection.Rnw:700-701
###################################################
  m <- selection(ys~xs, yo ~xo)


###################################################
### code chunk number 5: selection.Rnw:719-730
###################################################
par(mar=c(3,3,0,0) + 0.1,
    mgp=c(2,1,0))
pch <- c(1, 16)
plot(xo, yoX, pch=pch[1 + ys], cex=0.5, lwd=0.3)
abline(a=0, b=1, lty=1)
# True dependence
abline(a=coef(m)[3], b=coef(m)[4], lty=2)
# Heckman's model
cf <- coef(lm(yo ~ xo, subset=ys==1))
abline(a=cf[1], b=cf[2], lty=3)
# OLS


###################################################
### code chunk number 6: selection.Rnw:742-745
###################################################
yoX <- xs + eps[,2]
yo <- yoX*(ys > 0)
summary(selection(ys ~ xs, yo ~ xs))


###################################################
### code chunk number 7: selection.Rnw:764-775
###################################################
par(mar=c(3,3,0,0) + 0.1,
    mgp=c(2,1,0))
pch <- c(1, 16)
plot(xs, yoX, pch=pch[1 + ys], cex=0.5, lwd=0.3)
abline(a=0, b=1, lty=1)
# True dependence
abline(a=coef(m)[3], b=coef(m)[4], lty=2)
# Heckman's model
cf <- coef(lm(yo ~ xs, subset=ys==1))
abline(a=cf[1], b=cf[2], lty=3)
# OLS


###################################################
### code chunk number 8: selection.Rnw:801-806
###################################################
  xs <- runif(500, -5, 5)
  ys <- xs + eps[,1] > 0
yoX <- xs + eps[,2]
  yo <- yoX*(ys > 0)
  summary( selection(ys ~ xs, yo ~ xs))


###################################################
### code chunk number 9: selection.Rnw:808-809
###################################################
  m <- selection(ys ~ xs, yo ~ xs)


###################################################
### code chunk number 10: selection.Rnw:822-833
###################################################
par(mar=c(3,3,0,0) + 0.1,
    mgp=c(2,1,0))
pch <- c(1, 16)
plot(xs, yoX, pch=pch[1 + ys], cex=0.5, lwd=0.3)
abline(a=0, b=1, lty=1)
# True dependence
abline(a=coef(m)[3], b=coef(m)[4], lty=2)
# Heckman's model
cf <- coef(lm(yo ~ xs, subset=ys==1))
abline(a=cf[1], b=cf[2], lty=3)
# OLS


###################################################
### code chunk number 11: selection.Rnw:848-859
###################################################
  set.seed(0)
  vc <- diag(3)
  vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
  vc[upper.tri(vc)] <- vc[lower.tri(vc)]
  eps <- rmvnorm(500, c(0,0,0), vc)
  xs <- runif(500)
  ys <- xs + eps[,1] > 0
  xo1 <- runif(500)
  yo1 <- xo1 + eps[,2]
  xo2 <- runif(500)
  yo2 <- xo2 + eps[,3]


###################################################
### code chunk number 12: selection.Rnw:874-875
###################################################
  summary(selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2)))


###################################################
### code chunk number 13: selection.Rnw:885-896
###################################################
set.seed(5)
eps <- rmvnorm(1000, rep(0, 3), vc)
eps <- eps^2 - 1
xs <- runif(1000, -1, 0)
ys <- xs + eps[,1] > 0
xo1 <- runif(1000)
yo1 <- xo1 + eps[,2]
xo2 <- runif(1000)
yo2 <- xo2 + eps[,3]

summary(selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), iterlim=20))


###################################################
### code chunk number 14: code:t5chi_woER
###################################################
set.seed(6)
xs <- runif(1000, -1, 1)
  ys <- xs + eps[,1] > 0
  yo1 <- xs + eps[,2]
  yo2 <- xs + eps[,3]
summary(tmp <- selection(ys~xs, list(yo1 ~ xs, yo2 ~ xs), iterlim=20))


###################################################
### code chunk number 15: selection.Rnw:928-989
###################################################
   EUlower <- function(alpha) {
      alpha[alpha >= 1] <- NA
      alpha <- sqrt(-alpha + 1)
      EUalpha <- ss^2 - 2*ss*alpha*dnorm(alpha/ss)/(1 - 2*pnorm(-alpha/ss))
      s1s^2/ss^2*EUalpha + s1^2 - s1s^2/ss^2 - 1
   }
   EUupper <- function(alpha) {
      alpha[alpha >= 1] <- 1
      alpha <- sqrt(-alpha + 1)
      EUalpha <- (ss*alpha*dnorm(alpha/ss) + ss^2*(1 - pnorm(alpha/ss)))/(1 - pnorm(alpha/ss))
      s2s^2/ss^2*EUalpha + s2^2 - s2s^2/ss^2 - 1
   }
   Nlower <- function(alpha) {
      alpha <- -alpha
      -Ns1s*dnorm(alpha)/pnorm(alpha)
   }
   Nupper <- function(alpha) {
      alpha <- -alpha
      Ns2s*dnorm(-alpha)/pnorm(-alpha)
   }
   ss <- sqrt(vc[1,1])
   s1 <- sqrt(vc[2,2])
   s2 <- sqrt(vc[3,3])
   s1s <- vc[1,2]
   s2s <- vc[1,3]
   Ns1s <- coef(tmp)["sigma1"]*coef(tmp)["rho1"]
   Ns2s <- coef(tmp)["sigma2"]*coef(tmp)["rho2"]
hatb1O <- coef(tmp)[c("XO1(Intercept)", "XO1xs")]
hatb2O <- coef(tmp)[c("XO2(Intercept)", "XO2xs")]
   es <- eps[,1]
   e1 <- eps[,2]
e2 <- eps[,3]
   ex <- seq(-5, 5, length=200)
ey <- cbind(EUlower(ex), EUupper(ex),
##            -31*Nupper(cbind(1, ex)%*%hatb1O), 5.5*Nlower(cbind(1,ex)%*%hatb2O),
            -s1s*dnorm(-ex)/pnorm(-ex), 
            s2s*dnorm(ex)/pnorm(ex)
            )
   mcy <- matrix(0, length(ex), 2)
   for(i in seq(length=nrow(mcy))) {
      mcy[i,1] <- mean(e1[es < -ex[i]])
      mcy[i,2] <- mean(e2[es > -ex[i]])
   }
par(cex=0.8, mar=c(3,3,0,0) + 0.1, mgp=c(2,1,0))
matplot(ex, ey, type="l", lty=c(1,1,2,2,3,3), col=1,
        xlab=expression(x^S), ylab="",
        ylim=c(-1.5,2.5))
abline(v=-1,lty=3)
abline(v=1, lty=3)
#   matpoints(ex, mcy, pch=c(1,2), cex=0.5, col=1)
   axis(4)
text(-3.5, 0.8,
     expression(paste("correct  ")*E*group("[", epsilon^{O2}*group("|", epsilon^S > -bold(beta)^S*minute*bold(x)^S, ""), "]")))
text(-3.6, -0.4,
     expression(paste("correct  ")*E*group("[", epsilon^{O1}*group("|", epsilon^S < -bold(beta)^S*minute*bold(x)^S, ""), "]")))
text(-2.5, 1.5,
     expression(paste("assumed  ")*E*group("[", epsilon^{O2}*group("|",
         epsilon^S > -bold(beta)^S*minute*bold(x)^S, ""), "]")), pos = 4 )
text(1, -1.4,
     expression(paste("assumed  ")*E*group("[", epsilon^{O1}*group("|",
         epsilon^S < -bold(beta)^S*minute*bold(x)^S, ""), "]")), pos = 2 )


###################################################
### code chunk number 16: selection.Rnw:1011-1013
###################################################
  coef(summary(lm(yo1~xs, subset=ys==0)))
  coef(summary(lm(yo2~xs, subset=ys==1)))


###################################################
### code chunk number 17: greene22.8start
###################################################
data( "Mroz87" )
Mroz87$kids <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )


###################################################
### code chunk number 18: greene22.8TwoStep
###################################################
greeneTS <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city,
   data = Mroz87, method = "2step" )


###################################################
### code chunk number 19: greene22.8ML
###################################################
greeneML <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, data = Mroz87,
   maxMethod = "BHHH", iterlim = 500 )


###################################################
### code chunk number 20: cameron16.6start
###################################################
data( "RandHIE" )
subsample <- RandHIE$year == 2 & !is.na( RandHIE$educdec )
selectEq <- binexp ~ logc + idp + lpi + fmde + physlm + disea +
   hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
   child + fchild + black
outcomeEq <- lnmeddol ~ logc + idp + lpi + fmde + physlm + disea +
   hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
   child + fchild + black


###################################################
### code chunk number 21: cameron16.6TwoStep
###################################################
rhieTS <- selection( selectEq, outcomeEq, data = RandHIE[ subsample, ],
                          method = "2step" )


###################################################
### code chunk number 22: cameron16.6ML
###################################################
rhieML <- selection( selectEq, outcomeEq, data = RandHIE[ subsample, ] )


###################################################
### code chunk number 23: Greene22.8NoConvergence
###################################################
greeneStart <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city,
   data = Mroz87, start = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.9),
   control = list( qrtol = 1e-14 ) )

cat( greeneStart$message )

coef( summary( greeneStart ) )[ "rho", ]


###################################################
### code chunk number 24: Greene22.8SANN
###################################################
set.seed(0)
greeneSANN <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city,
   data = Mroz87, start = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.9),
   maxMethod="SANN", parscale = 0.001 )

greeneStartSANN <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city,
   data = Mroz87, start = coef( greeneSANN ) )

cat( greeneStartSANN$message )


###################################################
### code chunk number 25: selection.Rnw:1213-1215
###################################################
logLik( greeneML )
logLik( greeneStartSANN )


###################################################
### code chunk number 26: tobit_tobit2
###################################################
set.seed(0)
x <- runif(1000)
y <- x + rnorm(1000)
ys <- y > 0
tobitML <- selection( ys~x, y~x, control = list( qrtol = 1e-14 ) )
cat( tobitML$message )
coef( summary( tobitML ) )[ "rho", ]


###################################################
### code chunk number 27: tobit_tobit2_summary
###################################################
tobitTS <- selection(ys~x, y~x, method="2step")
coef( summary( tobitTS ) )[ "rho", ]


