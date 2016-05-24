### R code from vignette source 'pcse.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, digits = 3, useFancyQuotes = FALSE)
library("pcse")


###################################################
### code chunk number 2: pcse.Rnw:323-325
###################################################
library("pcse")
data("agl")


###################################################
### code chunk number 3: pcse.Rnw:328-330
###################################################
agl.lm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central +
  leftc + inter + as.factor(year), data = agl)


###################################################
### code chunk number 4: pcse.Rnw:344-345
###################################################
summary(agl.lm)


###################################################
### code chunk number 5: pcse.Rnw:348-349
###################################################
agl.pcse <- pcse(agl.lm, groupN = agl$country, groupT = agl$year)


###################################################
### code chunk number 6: pcse.Rnw:353-354
###################################################
summary(agl.pcse)


###################################################
### code chunk number 7: pcse.Rnw:365-371
###################################################
data("aglUn")
aglUn.lm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central +
  leftc + inter + as.factor(year), data = aglUn)
aglUn.pcse1 <- pcse(aglUn.lm, groupN = aglUn$country,
  groupT = aglUn$year, pairwise = TRUE)
summary(aglUn.pcse1)


###################################################
### code chunk number 8: pcse.Rnw:380-382
###################################################
aglUn.pcse2 <- pcse(aglUn.lm, groupN = aglUn$country, 
  groupT = aglUn$year, pairwise = FALSE)


###################################################
### code chunk number 9: pcse.Rnw:391-392
###################################################
summary(aglUn.pcse2)


