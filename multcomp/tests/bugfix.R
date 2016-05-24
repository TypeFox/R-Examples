
library("multcomp")
set.seed(290875)

### mcp didn't accept objects of class `contrMat'
### spotted by Yves Brostaux <brostaux.y@fsagx.ac.be>
amod <- aov(response ~ trt, data = cholesterol)
cht1 <- glht(amod, linfct = mcp(trt = "Tukey"))
K <- contrMat(table(cholesterol$trt), type = "Tukey")
cht2 <- glht(amod, linfct = mcp(trt = K))
stopifnot(all.equal(coef(cht1), coef(cht2)))

### several inconsistencies spotted by 
### Rich Heiberger <rmh@temple.edu> 2006-11-28

### need to be identical
stopifnot(identical(cht1, print(cht1)))

### was: error
summary(cht1)$test


### NAs in coefficients
tmp.data <- data.frame(EE=gl(2, 1, 24, letters[1:2]),
                FF=gl(3, 2, 24, LETTERS[3:5]),
                GG=gl(4, 6, 24, letters[6:9]))
tmp.data$x <- rep(12, 24)
tmp.data$y <- rep(7, 24)
tmp.data$z <- c(9, 14, 3, 4, 15, 1, 11, 13, 24, 10, 22, 18,
                20, 21, 6, 7, 16, 2, 19, 12, 17, 8, 23, 5)
tmp.data$w <- c(15, 9, 18, 21, 17, 11, 23, 12, 1, 10, 2, 14, 24, 7,
                13, 4, 5, 19, 16, 20, 3, 8, 22, 6)

tmp.aov <- aov(z ~ EE+FF*GG + x*y +x*EE + y*FF, data=tmp.data)

try(glht(tmp.aov, linfct=mcp(EE="Tukey")))
try(glht(tmp.aov, linfct=mcp(FF="Tukey")))
glht(tmp.aov, linfct=mcp(GG="Tukey"))

### covariate interactions: fire a warning
tmp.aov <- aov(z ~ w*GG , data=tmp.data)
glht(tmp.aov, linfct = mcp(GG = "Tukey"))

### stop with informative error message
amod <- aov(breaks ~ tension + Error(wool), data = warpbreaks)
try(glht(amod, linfct = mcp(tension = "Tukey")))

### print error, spotted by Rich
amod <- aov(breaks ~ wool * tension, data = warpbreaks)
wht <- glht(amod, linfct = mcp(tension = "Tukey"))
tmp <- confint(wht, calpha=2)
print(tmp)

### coef. and vcov. didn't pass through
### bug report by John Deke <jdeke73@gmail.com>
lmod <- lm(Fertility ~ ., data = swiss) 
my.model <- list(coef(lmod),vcov(lmod)) 
coef2 <- function(model) return(model[[1]]) 
vcov2 <- function(model) return(model[[2]]) 
a <- glht(model = my.model, linfct = c("Agriculture=0","Catholic=0"),
          coef. = coef2, vcov. = vcov2, df = 100) 
b <- glht(model = lmod, linfct = c("Agriculture=0","Catholic=0"), 
          df = 100)
stopifnot(all.equal(coef(a), coef(b)))

### checks in mcp (spotted by Rich)
amod <- aov(breaks ~ tension, data = warpbreaks)
try(glht(amod, linfct = mcp(group = "Tukey")))
tmp <- warpbreaks
class(tmp$tension) <- "numeric"
amod <- aov(breaks ~ tension, data = tmp)
try(glht(amod, linfct = mcp(tension = "Tukey")))

### symbolic description and interactions
### spotted by Antonio Fabio Di Narzo <antonio.dinarzo@unibo.it>
dat <- data.frame(y = rnorm(6), x = seq_len(6), f = gl(2, 3))
lf <- glht(lm(y ~ x * f, data = dat), 'x + x:f2 = 0')$linfct
stopifnot(all.equal(max(abs(lf - c(0, 1, 0, 1))), 0))
lf <- glht(lm(y ~ x * f, data = dat), 'x + 2.5 * x:f2 = 0')$linfct
stopifnot(all.equal(max(abs(lf - c(0, 1, 0, 2.5))), 0))

### example from Bretz 2001 JSCS

`tmp` <-
structure(list(gr = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), .Label = c("1", 
"2", "3"), class = "factor"), age = c(39L, 40L, 41L, 41L, 45L, 
49L, 52L, 47L, 61L, 65L, 58L, 59L, 29L, 29L, 33L, 32L, 31L, 29L, 
29L, 30L, 21L, 28L, 23L, 35L, 38L, 38L, 43L, 39L, 38L, 42L, 43L, 
43L, 37L, 50L, 50L, 45L, 48L, 51L, 46L, 58L, 27L, 25L, 24L, 32L, 
23L, 25L, 32L, 18L, 19L, 26L, 33L, 27L, 33L, 25L, 42L, 35L, 35L, 
41L, 38L, 41L, 36L, 36L, 41L, 41L, 37L, 42L, 39L, 41L, 43L, 41L, 
48L, 47L, 53L, 49L, 54L, 48L, 49L, 47L, 52L, 58L, 62L, 65L, 62L, 
59L), y = c(4.62, 5.29, 5.52, 3.71, 4.02, 5.09, 2.7, 4.31, 2.7, 
3.03, 2.73, 3.67, 5.21, 5.17, 4.88, 4.5, 4.47, 5.12, 4.51, 4.85, 
5.22, 4.62, 5.07, 3.64, 3.64, 5.09, 4.61, 4.73, 4.58, 5.12, 3.89, 
4.62, 4.3, 2.7, 3.5, 5.06, 4.06, 4.51, 4.66, 2.88, 5.29, 3.67, 
5.82, 4.77, 5.71, 4.47, 4.55, 4.61, 5.86, 5.2, 4.44, 5.52, 4.97, 
4.99, 4.89, 4.09, 4.24, 3.88, 4.85, 4.79, 4.36, 4.02, 3.77, 4.22, 
4.94, 4.04, 4.51, 4.06, 4.02, 4.99, 3.86, 4.68, 4.74, 3.76, 3.98, 
5, 3.31, 3.11, 4.76, 3.95, 4.6, 4.83, 3.18, 3.03)), .Names = c("gr", 
"age", "y"), row.names = c(NA, -84L), class = "data.frame")

amod <- aov(y ~ gr + age, data = tmp)
glht(amod, linfct = mcp(gr = "Tukey"))

### better error message
### suggested by Rich
amod <- aov(breaks ~ tension, data = warpbreaks)
try(glht(amod, linfct = mcp(tension = "Warp")))

### cld did not find a terms component
### spotted by Peter B. Mandeville <mandevip@hotmail.com>
if (require("nlme")) {
    data("Orthodont")
    fm1 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
    hsd1 <- glht(fm1, linfct = mcp(Sex = "Tukey"))
    cld(hsd1)
}

### spotted by <chris.chizinski@gmail.com>
### example code by Achim Zeileis <Achim.Zeileis@wu.ac.at>
## various models with and without intercept
m1a <- lm(breaks ~ tension, data = warpbreaks)
m1b <- lm(breaks ~ 0 + tension, data = warpbreaks)
m2a <- lm(breaks ~ wool + tension, data = warpbreaks)
m2b <- lm(breaks ~ 0 + wool + tension, data = warpbreaks)

## these two are equivalent: one factor with/without intercept
stopifnot(all.equal(
coef(glht(m1a, linfct = mcp(tension = "Tukey"))),
coef(glht(m1b, linfct = mcp(tension = "Tukey")))))

## these two should be equivalent: two factors with/without intercept
## but the latter fails
stopifnot(all.equal(
coef(glht(m2a, linfct = mcp(tension = "Tukey"))),
coef(glht(m2b, linfct = mcp(tension = "Tukey")))))

library("MASS")
xdf <- data.frame(y = gl(3, 10, ordered = TRUE), grp = sample(gl(3, 10)))
glht(polr(y ~ grp, data = xdf), mcp(grp = "Dunnett"))

### interactions of two factors
dat <- expand.grid(f = gl(2, 3), f2 = gl(3, 2))
dat$y <- rnorm(nrow(dat))
lf <- glht(lm(y ~ f : f2 - 1, data = dat), 'f1:f21 - f2:f22 = 0')$linfct
stopifnot(all.equal(max(abs(lf - c(1, 0, 0, -1, 0, 0))), 0))

### plotting one-sided confidence intervals
amod <- aov(breaks ~ wool + tension, data = warpbreaks)
wht <- glht(amod, linfct = mcp(tension = "Tukey"), alternative="greater")
plot(wht, xlim=c(-30, 30), main="right side was missing")
wht <- glht(amod, linfct = mcp(tension = "Tukey"), alternative="less")
plot(wht, xlim=c(-40, 20), main="left side was missing")

### reported by Christian Ritz
summary(glht(parm(1:4,matrix(c(1,0.97,0.89,0.74,
                               0.97,1,0.97,0.89,
                               0.89,0.97,1,0.97,
                               0.74,0.89,0.97,1), 4, 4))))

                               

library("multcomp")
set.seed(2343)

X <- data.frame(X1 = rep(c(1,0),c(20,30)),
  X2 = rep(rep(c(1,0),3),c(rep(10,4),0,10)),
  X3 = rep(rep(c(1,0),5),each=5))
Y <- rnorm(50,4 + 4*X[,1] + 4*X[,2] + X[,3] + .5*X[,1]*X[,3] + .4*X[,2]*X[,3],.25)

model <- lm(Y ~ (X1 + X2) * X3,data=X)
coef(model)

my.contrasts<- c(
  "X1 - X2 + .5*X1:X3 - .5*X2:X3 = 0",  # previously wrong answer (actually got X1 + X2 + 0.5* X1)
  "X1 + .5*X1:X3 - X2 - .5*X2:X3 = 0",  # previously wrong answer
  "X1 + .5*X1:X3 - .5*X2:X3 - X2 = 0")  # right answer

(contrast.result <- glht(model,lin = my.contrasts))

# right calculation
(ok <- sum(coef(model) * c(0,1,-1,0,.5,-.5)))

stopifnot(all.equal(as.numeric(coef(contrast.result)), rep(sum(coef(model) * c(0,1,-1,0,.5,-.5)),3)))


# actual calculation - note that -1 has changed to 1
#sum(coef(model) * c(0, 1, 1, 0, .5, -.5))


(mc <- multcomp:::chrlinfct2matrix(my.contrasts, names(coef(model))))

stopifnot(all.equal(as.numeric(mc$K[1,c('(Intercept)', 'X1','X2', 'X3','X1:X3','X2:X3')]),c( 0,1,-1 ,0, 0.5,-0.5)))
stopifnot(all.equal(as.numeric(mc$K[2,c('(Intercept)', 'X1','X2', 'X3','X1:X3','X2:X3')]),c( 0,1,-1 ,0, 0.5,-0.5)))
stopifnot(all.equal(as.numeric(mc$K[3,c('(Intercept)', 'X1','X2', 'X3','X1:X3','X2:X3')]),c( 0,1,-1 ,0, 0.5,-0.5)))

