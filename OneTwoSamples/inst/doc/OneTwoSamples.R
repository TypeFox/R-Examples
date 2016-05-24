### R code from vignette source 'OneTwoSamples.Rnw'

###################################################
### code chunk number 1: label
###################################################
library('OneTwoSamples')


###################################################
### code chunk number 2: label
###################################################
?data_outline
?t.test()


###################################################
### code chunk number 3: label
###################################################
## generate samples x and y
x = women$height; x
y = women$weight; y


###################################################
### code chunk number 4: label
###################################################
## operate on one sample
## one_two_sample(x) is equivalent to one_sample(x)
one_two_sample(x)


###################################################
### code chunk number 5: label
###################################################
## one_two_sample(y) is equivalent to one_sample(y)
one_two_sample(y)


###################################################
### code chunk number 6: hist_x
###################################################
x = women$height
## Histograms with density estimation curve and normal density curve
w<-seq(min(x),max(x),length.out = 51)
Vector = c(density(x)$y, dnorm(w, mean(x), sd(x)))
ylim = c(min(Vector), max(Vector))

hist(x, freq = FALSE, ylim = ylim, main = paste("Histogram of x"), xlab = "x")
lines(density(x),col="blue",lty = 1)
lines(w, dnorm(w, mean(x), sd(x)), col="red",lty = 2)
leg.txt = c("Density estimation curve","Normal density curve")
legend("topleft",legend = leg.txt,lty = 1:2,col = c('blue','red'))


###################################################
### code chunk number 7: hist_y
###################################################
y = women$weight
## Histograms with density estimation curve and normal density curve
w<-seq(min(y),max(y),length.out = 51)
Vector = c(density(y)$y, dnorm(w, mean(y), sd(y)))
ylim = c(min(Vector), max(Vector))

hist(y, freq = FALSE, ylim = ylim, main = paste("Histogram of y"), xlab = "y")
lines(density(y),col="blue",lty = 1)
lines(w, dnorm(w, mean(y), sd(y)), col="red",lty = 2)
leg.txt = c("Density estimation curve","Normal density curve")
legend("topleft",legend = leg.txt,lty = 1:2,col = c('blue','red'))


###################################################
### code chunk number 8: ecdf_x
###################################################
## Empirical cumulative distribution function (ECDF) vs normal cdf
plot(ecdf(x),verticals = TRUE, do.p = FALSE, main = "ecdf(x)", xlab = "x", ylab = "Fn(x)")
w<-seq(min(x),max(x),length.out = 51)
lines(w, pnorm(w, mean(x), sd(x)), col="red")


###################################################
### code chunk number 9: ecdf_y
###################################################
## Empirical cumulative distribution function (ECDF) vs normal cdf
plot(ecdf(y),verticals = TRUE, do.p = FALSE, main = "ecdf(y)", xlab = "y", ylab = "Fn(y)")
w<-seq(min(y),max(y),length.out = 51)
lines(w, pnorm(w, mean(y), sd(y)), col="red")


###################################################
### code chunk number 10: QQplot_x
###################################################
## QQ plot
qqnorm(x); qqline(x)


###################################################
### code chunk number 11: QQplot_y
###################################################
## QQ plot
qqnorm(y); qqline(y)


###################################################
### code chunk number 12: label
###################################################
## operate on two samples
one_two_sample(x, y)


