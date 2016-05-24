### R code from vignette source 'mhurdle.rnw'

###################################################
### code chunk number 1: mhurdle.rnw:105-106
###################################################
options(prompt= "R> ", useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mhurdle.rnw:1618-1620
###################################################
library("Formula")
f <- Formula(y ~ x11 + x12  | x21 + x22 | x31 + x32)


###################################################
### code chunk number 3: mhurdle.rnw:1727-1728
###################################################
library("mhurdle")


###################################################
### code chunk number 4: mhurdle.rnw:1756-1760
###################################################
data("Comics", package = "mhurdle")
head(Comics, 3)
mean(Comics$comics == 0)
max(Comics$comics)


###################################################
### code chunk number 5: mhurdle.rnw:1827-1834
###################################################
Comics$incu <- with(Comics, income / cu)
Comics$incum <- with(Comics, incu / mean(incu))


m010 <- mhurdle(comics ~ 0 | log(incum) + I(log(incum)^2) +
                I(log(incum)^3) + age  + gender + educ +
                size| 0, data = Comics, dist = "n", method = 'bfgs')


###################################################
### code chunk number 6: mhurdle.rnw:1847-1850
###################################################
m110d <- mhurdle(comics ~ gender + educ + age |  log(incum) +
                 I(log(incum)^2) + I(log(incum)^3) + size | 0,
                 data = Comics, corr = "d", dist = "n", method = 'bfgs')


###################################################
### code chunk number 7: mhurdle.rnw:1857-1858
###################################################
m110i <- update(m110d, corr = NULL)


###################################################
### code chunk number 8: mhurdle.rnw:1867-1868
###################################################
m100d <- update(m110d, dist = "ln")


###################################################
### code chunk number 9: mhurdle.rnw:1873-1874
###################################################
m100i <- update(m100d, corr = NULL)


###################################################
### code chunk number 10: mhurdle.rnw:1881-1884
###################################################
m111dii <- mhurdle(comics ~ gender + educ  |  log(incum) +
                   I(log(incum)^2) + I(log(incum)^3) + size | age,
                   data = Comics, corr = "dii", dist = "n", method = 'bfgs')


###################################################
### code chunk number 11: mhurdle.rnw:1893-1894
###################################################
summary(m111dii)


###################################################
### code chunk number 12: mhurdle.rnw:1910-1915
###################################################
coef(m111dii, "h2")
coef(m110d, "h1")
coef(m110d, "sd")
coef(summary(m111dii), "h3")
vcov(m111dii, "h3")


###################################################
### code chunk number 13: mhurdle.rnw:1922-1924
###################################################
logLik(m110d)
logLik(m110d, naive = TRUE)


###################################################
### code chunk number 14: mhurdle.rnw:1933-1934
###################################################
head(fitted(m110d))


###################################################
### code chunk number 15: mhurdle.rnw:1941-1949
###################################################
predict(m110d,
        newdata = data.frame(
            comics = c(0, 1, 2),
            gender = c("female", "female", "male"),
            age = c(20, 18, 32),
            educ = c(10, 20, 5),
            incum = c(4, 8, 2),
            size = c(2, 1, 3)))


###################################################
### code chunk number 16: mhurdle.rnw:1968-1969
###################################################
rsq(m110d, type = "coefdet")


###################################################
### code chunk number 17: mhurdle.rnw:1981-1982
###################################################
rsq(m110d, type = "lratio", adj = TRUE)


###################################################
### code chunk number 18: mhurdle.rnw:1994-1995
###################################################
vuongtest(m110d, m111dii)


###################################################
### code chunk number 19: mhurdle.rnw:2018-2020
###################################################
vuongtest(m100d, m100i, type = 'nested', hyp = TRUE)
vuongtest(m100d, m100i, type = 'nested', hyp = FALSE)


###################################################
### code chunk number 20: mhurdle.rnw:2034-2035
###################################################
coef(summary(m100d), "corr")


###################################################
### code chunk number 21: mhurdle.rnw:2047-2050
###################################################
m010bis <- mhurdle(comics ~ 0 | log(incum) + I(log(incum)^2) +
                   I(log(incum)^3)  + gender + educ + age +
                   empl+area| 0, data = Comics, dist = "n", method = 'bfgs')


###################################################
### code chunk number 22: mhurdle.rnw:2057-2058
###################################################
vuongtest(m010, m010bis, type="overlapping")


###################################################
### code chunk number 23: mhurdle.rnw:2065-2066
###################################################
vuongtest(m010, m010bis, type="non-nested")


###################################################
### code chunk number 24: mhurdle.rnw:2076-2077
###################################################
vuongtest(m010bis, m010, type="overlapping", hyp=TRUE)


