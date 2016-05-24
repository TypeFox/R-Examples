### R code from vignette source 'sclero_tutorial.Rnw'

###################################################
### code chunk number 1: initial_settings
###################################################
options(width=80, prompt = " ", continue = " ")


###################################################
### code chunk number 2: sclero_tutorial.Rnw:66-68 (eval = FALSE)
###################################################
## library(devtools)
## install_github("MikkoVihtakari/sclero", dependencies = TRUE)


###################################################
### code chunk number 3: sclero_tutorial.Rnw:145-146
###################################################
library(sclero)


###################################################
### code chunk number 4: sclero_tutorial.Rnw:151-153
###################################################
path <- file.path(system.file("extdata", package = "sclero"), "shellspots.zip")
dat <- read.ijdata(path, scale = 0.7812, unit = "um") 


###################################################
### code chunk number 5: sclero_tutorial.Rnw:158-159
###################################################
summary(dat)


###################################################
### code chunk number 6: sclero_tutorial.Rnw:165-166
###################################################
order.ijdata(dat, print.order = TRUE)


###################################################
### code chunk number 7: sclero_tutorial.Rnw:171-173
###################################################
dat2 <- order.ijdata(dat, gbs = c(1,3,6:14,4,5,2))
order.ijdata(dat2, print.order = TRUE)


###################################################
### code chunk number 8: sclero_tutorial.Rnw:179-180
###################################################
shell <- convert.ijdata(dat)


###################################################
### code chunk number 9: plot1
###################################################
plot(shell)


###################################################
### code chunk number 10: plot2
###################################################
aligned <- spot.dist(shell)
plot(aligned)


###################################################
### code chunk number 11: sclero_tutorial.Rnw:223-224
###################################################
aligned


###################################################
### code chunk number 12: sclero_tutorial.Rnw:229-231 (eval = FALSE)
###################################################
## aligned$output ## Results shown above
## aligned$det.dat ## Results not shown here to save space


###################################################
### code chunk number 13: sclero_tutorial.Rnw:284-286 (eval = FALSE)
###################################################
## data(shellspots)
## shell <- convert.ijdata(shellspots)


###################################################
### code chunk number 14: sclero_tutorial.Rnw:292-294
###################################################
path <- file.path(system.file("extdata", package = "sclero"))
shellsizes <- assign.size(shell, path = path)


###################################################
### code chunk number 15: sclero_tutorial.Rnw:298-299
###################################################
head(shellsizes$spot.area$spot.dat[[1]])


###################################################
### code chunk number 16: sclero_tutorial.Rnw:307-309
###################################################
spotsizes <- spot.dist(shellsizes)
head(spotsizes$output[[1]])


###################################################
### code chunk number 17: plotsize
###################################################
plot(spotsizes, spot.size = "actual")


###################################################
### code chunk number 18: sclero_tutorial.Rnw:328-330
###################################################
data(barium)
head(barium)


###################################################
### code chunk number 19: sclero_tutorial.Rnw:335-336
###################################################
data(shellsizes)


###################################################
### code chunk number 20: plotvalues
###################################################
shellvalues <- assign.value(shellsizes, barium, value.name = "Ba/Ca")
plot(shellvalues, spot.size = "actual", spot.type = "value", main.type = "none")


###################################################
### code chunk number 21: plotvaluesalign
###################################################
shellvalues.aligned <- spot.dist(shellvalues)
plot(shellvalues.aligned, spot.size = "actual", spot.type = "idvalue",
  spot.color = "darkgrey", highlight.gbs = c("WG_start", "WG_end"))


###################################################
### code chunk number 22: sclero_tutorial.Rnw:379-385
###################################################
file <- file.path(system.file("extdata", package = "sclero"), "multi_spotseq.zip")
dat <- read.ijdata(file, scale = 0.7812, unit = "um")
multispot.raw <- convert.ijdata(dat)
path <- file.path(system.file("extdata", package = "sclero"))
multispot.size <- assign.size(multispot.raw, path = path)
multispot <- spot.dist(multispot.size)


###################################################
### code chunk number 23: multispot
###################################################
plot(multispot, spot.size = "actual")


