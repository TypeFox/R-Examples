## Test provided by Eddy Delpierre 2007-01-11


#data1 <- read.table("c://stat//projects//R//drcurves//r-part//pkgfolder//drc//tests//test1.data1.txt", header = TRUE)
data1 <- read.table("test1.data1.txt", header = TRUE)
#data2 <- read.table("c://stat//projects//R//drcurves//r-part//pkgfolder//drc//tests//test1.w1.txt", header = TRUE)
data2 <- read.table("test1.w1.txt", header = TRUE)

library(drc)

#FIT1 <- drm(y~x, data = data1, fct = LL.4())
#summary(FIT1)

FIT2 <- drm(y~x, data = data1, fct = LL.4(method = "3"))
summary(FIT2)

FIT3 <- drm(y~x, data = data1, fct = LL.4(method = "3"), weights = data2[, 2])
summary(FIT3)
plot(FIT3)

FIT4a <- drm(y~x, data=data1, fct=LL.4(fixed=c(6, -5e-9, NA, NA), method = "3"))
summary(FIT4a)
plot(FIT4a)

FIT4 <- drm(y~x, data=data1, fct=LL.4(fixed=c(6, -5e-9, NA, NA), method = "3"), weights=data2[, 2])
summary(FIT4)
plot(FIT4)

FIT5 <- drm(y~x, data=data1, fct=LL.4(fixed=c(NA,1E-9,5E-8,1E+12)))
summary(FIT5)

FIT6 <- drm(y~x, data=data1, fct=LL.4(fixed=c(NA,1E-9,5E-8,1E+12), method = "2"))
summary(FIT6)
