### R code from vignette source 'Rstyle.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Rstyle.Rnw:24-25
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: Rstyle.Rnw:167-168
###################################################
dir.create("plots", showWarnings=F)


###################################################
### code chunk number 3: Roptions
###################################################
options(width=100, continue="+ ")
options(useFancyQuotes = FALSE) 
set.seed(12345)
pdf.options(onefile=F,family="Times",pointsize=12)


###################################################
### code chunk number 4: Rstyle.Rnw:941-948
###################################################
x <- runif(1000, min = 0, max = 100)
xf <- cut(x, breaks = c(-1, 20, 50, 80, 101), labels = c("cold", "luke", "warm", "hot"))
xfdummies <- contrasts(xf, contrasts = FALSE )[xf,]
colnames(xfdummies) <-  paste("xf", c("cold", "luke", "warm", "hot"), sep="")
rownames(xfdummies) <- names(x)
dat <- data.frame(x, xf, xfdummies)
head(dat)


###################################################
### code chunk number 5: Rstyle.Rnw:956-967 (eval = FALSE)
###################################################
## set.seed(12345)
## x1 <- rnorm(200, m = 300, s = 140)
## x2 <- rnorm(200, m = 80, s = 30)
## y <- 3 + 0.2 * x1 + 0.4 * x2 + rnorm(200, s=400)
## dat <- data.frame(x1, x2, y); rm(x1,x2,y)
## m1 <- lm (y ~ x1 + x2, data = dat)
## m1summary <- summary(m1)
## m1se <- m1summary$sigma
## m1rsq <- m1summary$r.squared
## m1coef <- m1summary$coef
## m1aic <- AIC(m1)


###################################################
### code chunk number 6: ps10 (eval = FALSE)
###################################################
## library(rockchalk)
## dat$y2 = with(dat, 3 + 0.02 * x1 + 0.05 * x2 + 2.65 * x1 *x2 + rnorm(200, s=4000))
## par(mfcol=c(1,2))
## m1 <- lm(y2 ~ x1 + x2, data = dat)
## m1i <- lm(y2 ~ x1 * x2, data = dat)
## m1ps <- plotSlopes(m1, plotx = "x1", modx = "x2")
## m1ips <- plotSlopes(m1i, plotx = "x1", modx = "x2")
## m1imc <- meanCenter(m1i)
## m1irc <- residualCenter(m1i)	


