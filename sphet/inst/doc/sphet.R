### R code from vignette source 'sphet.Rnw'

###################################################
### code chunk number 1: optionsset
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: readdata
###################################################
library("sphet")
library(spdep)
data("boston", package = "spdep")


###################################################
### code chunk number 3: defformula
###################################################
listw <- nb2listw(boston.soi)


###################################################
### code chunk number 4: generatecoordinates
###################################################
set.seed(1234)
X <- runif(100, 0, 70)
Y <- runif(100, -30, 20)
coord1 <- cbind(seq(1, 100), X, Y)


###################################################
### code chunk number 5: distance1
###################################################
thm1 <- distance(coord = coord1, region.id = NULL, output = FALSE, type = "inverse", measure = "euclidean")
print(thm1[1:15, ])


###################################################
### code chunk number 6: distance2
###################################################
thm2 <- distance(coord1, region.id = NULL, output = FALSE, type = "NN", nn = 6)
print(thm2[1:15, ])


###################################################
### code chunk number 7: distance3
###################################################
thm3 <- distance(coord1, region.id = NULL, output = TRUE, type = "distance", cutoff = 1, measure = "gcircle", shape.name = "shapefile", region.id.name = "id1", firstline = TRUE, file.name = "dist_100.GWT")
class(thm3)


###################################################
### code chunk number 8: read
###################################################
id1 <- seq(1, nrow(boston.utm))
tmp <- distance(boston.utm, region.id = id1, output = TRUE, type = "NN", nn = 10, shape.name = "shapefile", region.id.name = "id1", firstline = TRUE, file.name = "boston_nn_10.GWT")
coldist <- read.gwt2dist(file = "boston_nn_10.GWT", region.id = id1, skip = 1)


###################################################
### code chunk number 9: read2
###################################################
class(coldist)


###################################################
### code chunk number 10: printsummarydist
###################################################
summary(coldist)


###################################################
### code chunk number 11: stslshac
###################################################
res <- stslshac(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw, distance = coldist, type = "Triangular", HAC = TRUE) 
summary(res)


###################################################
### code chunk number 12: stslshac1
###################################################
res <- stslshac(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw, distance = coldist, HAC = FALSE)
summary(res)


###################################################
### code chunk number 13: stsls
###################################################
res <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw = listw)
summary(res)


###################################################
### code chunk number 14: stsls2
###################################################
res <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw = listw, robust = TRUE)
summary(res)


###################################################
### code chunk number 15: stslshac2
###################################################
fix <- max(unlist(attributes(coldist)$GeoDa$dist))
res <- stslshac(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw, distance = coldist, type = "Parzen", bandwidth = fix)
summary(res)


###################################################
### code chunk number 16: gstslshet1
###################################################
res <- gstslshet(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw = listw, initial.value = 0.2)
summary(res)


###################################################
### code chunk number 17: gstslshet2
###################################################
res <- gstslshet(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c, listw = listw, initial.value = 0.2, inverse = FALSE, eps = 1e-18, sarar = FALSE )
summary(res)


