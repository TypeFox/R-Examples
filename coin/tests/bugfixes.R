### Regression tests for fixed bugs

set.seed(290875)
library("coin")
isequal <- coin:::isequal
GE <- coin:::GE
options(useFancyQuotes = FALSE)

### I() returns objects of class AsIs which caused an error in `trafo'
df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = gl(2, 50))
independence_test(I(x1 / x2) ~ x3, data = df)
independence_test(I(x1 < 0) ~ x3, data = df)

### expectation was wrong when varonly = TRUE in case both
### xtrafo and ytrafo were multivariate
if (require("multcomp")) {
    df <- data.frame(x = runif(30), y = runif(30), z = gl(3, 10))
    a <- independence_test(x + y ~ z, data = df,
         distribution = approximate(B = 19999),
         xtrafo = function(data) trafo(data, factor_trafo = function(x)
             model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))))
    b <- independence_test(x + y ~ z, data = df,
         xtrafo = function(data) trafo(data, factor_trafo = function(x)
             model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))))
    isequal(expectation(a), expectation(b))
}


### `statistic' for linear and standardized statistics was wrong in case of
### scores
data("jobsatisfaction")
stopifnot(unique(dim(statistic(lbl_test(jobsatisfaction), "linear"))) == 1)
stopifnot(unique(dim(statistic(lbl_test(jobsatisfaction), "standardized"))) == 1)


### support() failed in most cases
df <- data.frame(x = runif(20), y = runif(20), z = gl(2, 10))
support(independence_test(x ~ z, data = df))
support(independence_test(x ~ z, data = df, teststat = "quad"))
ite <- independence_test(I(round(x, 1)) ~ z, data = df, dist = exact())
ae <- support(ite)
de <- sapply(ae, function(x) dperm(ite, x))
sum(de)
ita <- independence_test(I(round(x, 1)) ~ z, data = df,
                         dist = approximate(B = 100000))
aa <- support(ita)
da <- sapply(aa, function(x) dperm(ita, x))
sum(da)
mean(round(ae, 10) %in% round(aa, 10))

plot(aa, da, type = "s", lty = 1)
lines(ae, de, type = "s", lty = 2)
itas <- independence_test(I(round(x, 1)) ~ z, data = df)
lines(ae[-1], diff(sapply(ae, function(x) pperm(itas, x))), lty = 3)
legend("topleft", lty = 1:3, legend = c("approx", "exact", "asympt"), bty = "n")

### check correct handling of multiple censoring indicators (in modeltools)
### was never wrong, just in case...
data("photocar", package = "coin")
i1 <- independence_test(Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                  data = photocar)
i2 <- independence_test(Surv(time, event) ~ group, data = photocar)
i3 <- independence_test(Surv(dmin, tumor) ~ group, data = photocar)

stopifnot(max(abs(statistic(i1, "standardized")[,1] -
                  statistic(i2, "stand"))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(statistic(i1, "standardized")[,2] -
                  statistic(i3, "stand"))) < sqrt(.Machine$double.eps))

### check new var_trafo argument
x <- rnorm(20)
y <- gl(2, 10)
a <- trafo(data.frame(x = x, y = y), numeric_trafo = normal_trafo)
b <- trafo(data.frame(x = x, y = y), var_trafo = list(x = normal_trafo))
stopifnot(isequal(a, b))

### check for multiple ordered factors
mydf <- data.frame(x = ordered(gl(4, 5)), y = ordered(gl(5, 4)),
                   z = rnorm(20))
it1 <- independence_test(x + y ~ z , data = mydf)
stopifnot(isequal(drop(statistic(it1, "linear")),
          c(statistic(independence_test(x ~ z , data = mydf), "linear"),
            statistic(independence_test(y ~ z , data = mydf), "linear"))))
it1 <- independence_test(x + z ~ y , data = mydf)
stopifnot(isequal(drop(statistic(it1, "linear")),
          c(statistic(independence_test(x ~ y , data = mydf), "linear"),
            statistic(independence_test(z ~ y , data = mydf), "linear"))))
it1 <- independence_test(z ~ x + y , data = mydf)
stopifnot(isequal(drop(statistic(it1, "linear")),
          c(statistic(independence_test(z ~ x , data = mydf), "linear"),
            statistic(independence_test(z ~ y , data = mydf), "linear"))))

### NA's and weights
mydf <- data.frame(x = 1:10, y = gl(2, 5), w = rep(2, 10))
s <- statistic(independence_test(x ~ y, data = mydf, weights = ~ w), "linear")
stopifnot(s == 30)
mydf$x[1] <- NA
s <- statistic(independence_test(x ~ y, data = mydf, weights = ~ w), "linear")
stopifnot(s == 28)

### two observations only
mydf <- data.frame(x = 1:10, y = factor(rep(c(1, 2), 5)))
independence_test(y ~ x, data = mydf, subset = c(1, 6))
independence_test(y ~ x, data = mydf, subset = c(1, 2))
try(independence_test(y ~ x, data = mydf, subset = 1))

### names of expectation and covariance
YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

it <- independence_test(length ~ site, data = YOY,
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
    teststat = "quad")
expectation(it)
covariance(it)

mydf <- data.frame(x = rnorm(10), y = rnorm(10), z = gl(2, 5))
it <- independence_test(x + y ~ z, data = mydf)
statistic(it, "linear")
expectation(it)
covariance(it)

### maxstat_trafo
n <- seq(from = 5, to = 100, by = 1)
for (i in n) {
   x <- round(rnorm(i) * 10, 1)
   xm <- maxstat_trafo(x)
   stopifnot(min(c(mean(xm[,1]), 1 - mean(xm[,ncol(xm)])) - 0.1) >
             -sqrt(.Machine$double.eps))
}

### formula evaluation in `parent.frame()', spotted by Z
foo <- function(x, y) independence_test(y ~ x)
a <- 1:10
b <- 1:10
foo(a, b)
x <- 1
y <- 1
foo(a, b)

### factors with only one level
dat <- data.frame(y = rnorm(100), x1 = runif(100), x2 = factor(rep(0, 100)))
try(independence_test(y ~ x1  + x2, data = dat))

### user specified g: names, MC
me <- as.table(matrix(c( 6,  8, 10,
               32, 47, 20), byrow = TRUE, nrow = 2,
    dimnames = list(group = c("In situ", "Control"),
                    genotype = c("AA", "AG", "GG"))))
medf <- as.data.frame(me)

add <- c(0, 1, 2)
dom <- c(0, 1, 1)
rez <- c(0, 0, 1)
g <- function(x) {
    x <- unlist(x)
    cbind(add[x], dom[x], rez[x])
}
it <- independence_test(group ~ genotype,
    data = medf, weights = ~ Freq, xtrafo = g)
statistic(it, "linear")

it <- independence_test(group ~ genotype,
    data = medf, weights = ~ Freq, xtrafo = g,
    distribution = approximate(B = 49999))
pvalue(it)

stopifnot(all.equal(statistic(independence_test(t(me), xtrafo = g), "linear"),
                    statistic(it, "linear")))

### alternative trafo for ordered variables didn't work
### spotted by Ludwig Hothorn <hothorn@biostat.uni-hannover.de>
tmp <- data.frame(x = ordered(sample(1:3, 20, replace = TRUE)), y = rnorm(20))
it1 <- independence_test(y ~ x, data = tmp, scores = list(x = c(1, 1, 2)))
g <- function(x) c(1, 1, 2)[unlist(x)]
it2 <- independence_test(y ~ x, data = tmp, xtrafo = g)
it3 <- independence_test(y ~ x, data = tmp,
    xtrafo = function(data) trafo(data, ordered_trafo = g))
stopifnot(all.equal(statistic(it1), statistic(it2)))
stopifnot(all.equal(statistic(it1), statistic(it3)))

### precision problems in SR algorithm, >= <= issue
### spotted by "Fay, Michael (NIH/NIAID) [E]" <mfay@niaid.nih.gov>
x <- c(1,2,3.1,4,5,6)
g <- factor(c(0,0,0,1,1,1))
it <- independence_test(x ~ g, distribution = exact())
stopifnot(pvalue(it) == 0.1)
itMC <- independence_test(x ~ g, distribution = approximate(99999))
ci <- attr(pvalue(itMC), "conf.int")
stopifnot(ci[1] < pvalue(it) && ci[2] > pvalue(it))

### any() not applied to logicals
### spotted by R-2.7.0 and Kaspar Daniel Hansen
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
### must not give a warning
wilcoxsign_test(x ~ y, alternative = "greater",
                distribution = exact())

### inconsistencies with confidence intervals
### spotted by Fritz Scholz <fscholz@u.washington.edu>
Route = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L,
2L, 2L, 2L), .Label = c("A", "B"), class = "factor")
Route.Time = c(5.8, 5.8, 5.9, 6, 6, 6, 6.3, 6.3, 6.4, 6.5, 6.5, 6.5,
6.8, 7.1, 7.3, 10.2)
Route2 <- factor(as.character(Route), levels = c("B", "A"))

wilcox_test(Route.Time~Route,conf.int=TRUE)
wilcox_test(Route.Time~Route2,conf.int=TRUE)
wilcox_test(Route.Time~Route,conf.int=TRUE, distr = exact())
wilcox_test(Route.Time~Route2,conf.int=TRUE, distr = exact())

### evaluate all formulae in xxx_tests parent.frame
### spotted by Matthieu Dubois <matthdub@gmail.com>
y <- as.data.frame(matrix(rnorm(200), ncol=2))
group <- gl(2, 50)
lapply(y, function(var) wilcox_test(var ~ group))

### make sure a parallel design with
### n = 2 isn't confused with a block design
### spotted by Fritz Scholz <fscholz@u.washington.edu>
dat <- data.frame(y = c(1, 2), g = gl(2, 1))
wt <- wilcox_test(y ~ g, data = dat, distribution = exact())
s <- support(wt)
stopifnot(all(dperm(wt, s) == c(0.5, 0.5)))

### dpqperm
wte <- wilcox_test(Route.Time ~ Route, distribution = exact())
wta <- wilcox_test(Route.Time ~ Route, distribution = approximate())
de <- dperm(wte, support(wte))
pe <- pperm(wte, support(wte))
stopifnot(max(abs(cumsum(de) - pe)) < sqrt(.Machine$double.eps))
da <- dperm(wta, support(wta))
pa <- pperm(wta, support(wta))
stopifnot(max(abs(cumsum(da) - pa)) < sqrt(.Machine$double.eps))
qperm(wte, seq(from = 0.1, to = 0.9, by = 0.1))
qperm(wta, seq(from = 0.1, to = 0.9, by = 0.1))

### split-up and confint (spotted by Fritz Scholz <fscholz@u.washington.edu>)
iseed <- 25
n <- m <- 4
conf.level <- 0.98
alternative <- "greater"
set.seed(iseed)
x <- round(rnorm(m), 2)
y <- round(rnorm(n, 2), 2)
score.factor <- factor(c(rep("Y",n),rep("X",m)),
                       levels = c("Y","X"))
dat.fr <- data.frame(scores=c(y,x), fac = score.factor)
it <- normal_test(scores ~ fac, data = dat.fr,
    conf.int = TRUE, conf.level = conf.level,
    alternative = alternative, dist = exact())
confint(it)

### discrete (spotted by Henric Winell <henric.winell@sorch.se>)
set.seed(1)
x <- gl(3, 5)
y1 <- rnorm(15)
y2 <- rnorm(15)
it <- independence_test(y1 + y2 ~ x, distribution = approximate(B = 5))
pvalue(it, "discrete") # didn't work

### error messages
### first group completely empty
x <- c(NA, NA, NA)
y <- c(2,4, 3)
z <- factor(c(1,1,1,2,2,2))
u <- c(x,y)
try(wilcox_test(u ~ z))

### missing values
x <- c(NA, NA, 1)
y <- c(2,  NA, NA)
z <- factor(c(1,1,1,2,2,2))
u <- c(x,y)
wilcox_test(u ~ z)

x <- c(NA, 1, 2, 8)
y <- c(2, 4, NA, 3)
wilcoxsign_test(x ~ y)

### no observations left after removing NAs
x <- c(NA, NA)
y <- c(1, 1)
try(wilcoxsign_test(x ~ y))

### problem in coin 1.0-12 fixed
x <- c(1, 2, 3)
y <- c(0, 1, 0)
wilcoxsign_test(x ~ y)

### as.integer(round(x)) is safer than as.integer(x)
 water_transfer <- data.frame(
     pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91,
            1.64, 0.73, 1.46, 1.15, 0.88, 0.90, 0.74, 1.21) * 100 - 72,
     age = factor(c(rep("At term", 10),
                    rep("12-26 Weeks", 5))))

### check this out
# water_transfer$pd[11]
# as.integer(water_transfer$pd[11])
# as.integer(round(water_transfer$pd[11]))

p1 <- pvalue(independence_test(pd ~ age, data = water_transfer,
                               alternative = "less",
                               distribution = exact(algorithm = "shift")))
p2 <- pvalue(independence_test(pd ~ age, data = water_transfer,
                               alternative = "less",
                               distribution = exact(algorithm = "split")))
stopifnot(isequal(p1, p2))

### monotonicity wasn't enforced for "step-down"
set.seed(290875)

gr <- gl(2, 50)
x1 <- rnorm(100) + (as.numeric(gr) - 1) * 0.5
x2 <- rnorm(100) - (as.numeric(gr) - 1) * 0.5

it <- independence_test(x1 + x2 ~ gr, distribution = approximate(B = 1000))

psd <- pvalue(it, "step-down") # wasn't monotone
stopifnot(psd[1] == psd[2])

### single-step p-values were too small
df <- data.frame(y1 = c(6, 7, 8, 5, 4, 3, 1, 2),
                 y2 = c(1, 2, 5, 4, 7, 3, 8, 6),
                 y3 = c(5, 7, 8, 6, 2, 3, 1, 4),
                 y4 = c(4, 8, 7, 3, 6, 5, 1, 2),
                 x = gl(2, 4, labels = c("I", "II")))

set.seed(711109)
it <- independence_test(y1 + y2 + y3 + y4 ~ x, data = df,
                        alternative = "greater",
                        distribution = approximate(B = 20))

pss <- pvalue(it, "single-step")
psd <- pvalue(it, "step-down")
stopifnot(isequal(all(GE(pss, psd)), TRUE))

### fmaxstat_trafo 'drop'ed its dimensions
fmaxstat_trafo(gl(2, 2))

### fmaxstat_trafo didn't drop unused levels
a <- gl(4, 2)
a[5:6] <- NA
fmaxstat_trafo(a)

### wrong p-value in the one-sided case
median_test(y1 ~ x, data = df, alternative = "less")

### nullvalue never got printed
logrank_test(Surv(time, event) ~ stadium, data = ocarcinoma)

### partial matching of the scores argument didn't work
chisq_test(as.table(jobsatisfaction[, , "Female"]), sco = list(Income = 1:4))

### partial matching of the scores argument didn't work
tab <- as.table(matrix(c(1, 2, 3, 1, 2, 3, 1, 2, 3), ncol = 3))
mh_test(tab, sco = list(response = 1:3))

### erroneously claimed to be linear-by-linear tests
chisq_test(as.table(jobsatisfaction[, , "Female"]), scores = list(Income = 1:4))
cmh_test(as.table(jobsatisfaction[, , "Female"]), scores = list(Income = 1:4))

### multivariate trafos with blocks failed for ordered factors and survival data
z <- gl(2, 10)
trafo(data.frame(gl(4, 5, ordered = TRUE)),
      ordered_trafo = function(x) cbind(c(1, 2, 3, 4)[x], c(1, 2, 3, 5)[x]),
      block = z)

### wrong exact p-values for stratified tests when all blocks have length two
y <- 1:30
x <- factor(rep(1:2, 15))
id <- gl(15, 2)
it <- independence_test(y ~ x | id, distribution = "exact") # Wrong
wt <- wilcoxsign_test(y ~ x | id, distribution = "exact") # OK! p = 6.104e-5
stopifnot(isequal(pvalue(it), pvalue(wt)))

### chisq_test standardized test statistic was wrong
ct <- chisq_test(as.table(matrix(1:4, ncol = 2)))
stopifnot(isequal(statistic(ct), statistic(ct, "standardized")^2))

### standardized signed-rank statistic had wrong sign
y1 <- seq(1, 30, 2)
y2 <- seq(2, 30, 2)
dta <- stack(list(y1 = y1, y2 = y2))
dta$block <- factor(rep(seq_along(y1), 2))

wt0 <- wilcox.test(y1, y2, paired = TRUE, exact = FALSE, correct = FALSE,
                   alternative = "greater")
wt1 <- wilcoxsign_test(y1 ~ y2, alternative = "greater")
wt2 <- wilcoxsign_test(values ~ ind | block, data = dta,
                       alternative = "greater")
it <- independence_test(values ~ ind | block, data = dta,
                        alternative = "less")

stopifnot(isequal(wt0$statistic, statistic(wt1, "linear"))) # Was OK
stopifnot(isequal(wt0$statistic, statistic(wt2, "linear")))
stopifnot(isequal(statistic(wt1), statistic(wt2)))
stopifnot(isequal(statistic(it), statistic(wt2)))

### maxstat_test didn't handle scores properly
mt0 <- maxstat_test(counts ~ age, data = treepipit)

fage <- factor(treepipit$age) # -> wrong estimates
mt1 <- maxstat_test(counts ~ fage, data = treepipit,
                    scores = list(fage = 1:4))
stopifnot(isequal(mt0@estimates$estimate$cutpoint,
                  max(as.numeric(mt1@estimates$estimate$cutpoint))))

oage <- ordered(treepipit$age) # -> error: oage is not a factor
mt2 <- maxstat_test(counts ~ oage, data = treepipit,
                    scores = list(oage = 1:4))
stopifnot(isequal(mt0@estimates$estimate$cutpoint,
                  max(as.numeric(mt2@estimates$estimate$cutpoint))))

### one-sided Ansari-Bradley test reported the wrong direction
y <- c(rnorm(10, sd = 1),  rnorm(10, sd = 5)) # sigma_1 < sigma_2 ==>> "less"
x <- gl(2, 10)
alt <- "less"
at <- ansari_test(y ~ x, alternative = alt)
stopifnot(isequal(at@statistic@alternative, alt))

### objects of class "SymmetryProblem" didn't check validity
dta <- data.frame(y = rnorm(20), y2 = rnorm(20),
                  x = factor(rep(1:4, 5)), x2 = gl(2, 10),
                  b = gl(5, 4),
                  w = rep(2, 20))
try(friedman_test(y + y2 ~ x | b, data = dta))
try(friedman_test(y ~ x + x2 | b, data = dta))
try(friedman_test(y ~ x2 | b, data = dta)) # was ok

### friedman_test didn't warn on weights
friedman_test(y ~ x | b, data = dta,  weights = ~ w)

### chisq_test ignored xtrafo and ytrafo
chisq_test(as.table(jobsatisfaction[, , "Female"]),
           xtrafo = function(data)
               trafo(data, factor_trafo = function(x)
                   of_trafo(x, scores = 1:4)))
chisq_test(as.table(jobsatisfaction[, , "Female"]),
           ytrafo = function(data)
               trafo(data, factor_trafo = function(y)
                   of_trafo(y, scores = 1:4)))
chisq_test(as.table(jobsatisfaction[, , "Female"]),
           xtrafo = function(data)
               trafo(data, factor_trafo = function(x)
                   of_trafo(x, scores = 1:4)),
           ytrafo = function(data)
               trafo(data, factor_trafo = function(y)
                   of_trafo(y, scores = 1:4)))

### dperm could fail for *exact* tests with blocks
### (this was due to the support values not being ordered and distinct)
dta <- data.frame(
    y = c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30,
          0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29),
    x = gl(2, 9),
    b = factor(rep(seq_len(9), 2)))
it <- independence_test(y ~ x | b, data = dta,
                        distribution = exact(algorithm = "shift"),
                        alternative = "greater")
stopifnot(is.numeric(dperm(it, support(it))))
### see 'regtest_distribution.R' for more extensive checks

### qperm with p = 1 could be NA for *exact* tests using the shift algorithm
TeaTasting <-
    matrix(c(3, 1, 1, 3),
           nrow = 2,
           dimnames = list(Guess = c("Milk", "Tea"),
                           Truth = c("Milk", "Tea")))
it <- independence_test(as.table(TeaTasting),
                        distribution = exact(algorithm = "shift"))
stopifnot(!is.na(qperm(it, p = 1)))

### dperm for *exact* tests returned a zero-length vector for values not
### included in the support
dp_it <- dperm(it, 1.5)
stopifnot(length(dp_it > 0) && dp_it == 0)

### dperm for *approximate* tests returned non-zero probabilities for values not
### included in the support
exfoliative <- matrix(c( 0, 16,
                        15, 57),
                      nrow = 2, byrow = TRUE)
set.seed(290875)
it <- independence_test(as.table(exfoliative),
                        distribution = approximate(B = 10000),
                        teststat = "scalar")
stopifnot(isequal(round(dperm(it, -statistic(it)), 4), 0.0000)) # 0.0747

### 'is_2sample' didn't drop unused levels, causing trouble with, e.g., subset
dta <- data.frame(y = rnorm(15), x = gl(3, 5), b = factor(rep(1:5, 3)))
subs <- dta$x %in% 1:2
wilcox_test(y ~ x | b, data = dta, subset = subs)

### problem with subsetting stratification variable
subs <- dta$x %in% 1:2 & dta$b %in% 1:4
wilcox_test(y ~ x | b, data = dta, subset = subs)

### 'dperm' returned non-sense for asymptotic max-type tests
y1 <- rnorm(10)
y2 <- rnorm(10)
x <- gl(2, 5)
it <- independence_test(y1 + y2 ~ x)
stopifnot(is.na(dperm(it, statistic(it))))

### 'p-'/'qperm' for asymptotic max-type test had issues with vector arguments
stopifnot(isequal(pperm(it, 1:3),
                  c(pperm(it, 1), pperm(it, 2), pperm(it, 3))))
stopifnot(isequal(qperm(it, c(0.9, 0.95, 0.99)),
                  c(qperm(it, 0.9), qperm(it, 0.95), qperm(it, 0.99))))

### blockwise permutations were only correct for factors ordered wrt their levels
set.seed(36)
.Call("R_blockperm", rep(1:4, 2),                       # was OK
      PACKAGE = "coin")

set.seed(36)
.Call("R_blockperm", rep(1:4, each = 2),                # was OK
      PACKAGE = "coin")

set.seed(36)
.Call("R_blockperm", c(1:4, 4:1),                       # was OK
      PACKAGE = "coin")

set.seed(36)
.Call("R_blockperm", c(4L, 1L, 2L, 2L, 4L, 3L, 1L, 3L), # wrong
      PACKAGE = "coin")

### could not distinguish censored and numeric responses
y1 <- rnorm(10)
y2 <- rnorm(10)
x <- gl(2, 5)
b <- factor(rep(1:5, 2))
try(oneway_test(Surv(y1) ~ x))
try(wilcox_test(Surv(y1) ~ x))
try(kruskal_test(Surv(y1) ~ x))
try(normal_test(Surv(y1) ~ x))
try(median_test(Surv(y1) ~ x))
try(savage_test(Surv(y1) ~ x))
try(kruskal_test(Surv(y1) ~ x))
try(taha_test(Surv(y1) ~ x))
try(klotz_test(Surv(y1) ~ x))
try(mood_test(Surv(y1) ~ x))
try(ansari_test(Surv(y1) ~ x))
try(fligner_test(Surv(y1) ~ x))
try(conover_test(Surv(y1) ~ x))
try(sign_test(Surv(y1) ~ x | b))
try(sign_test(Surv(y1) ~ Surv(y2)))
try(sign_test(Surv(y1) ~ y2))
try(sign_test(y1 ~ Surv(y2)))
try(wilcoxsign_test(Surv(y1) ~ x | b))
try(wilcoxsign_test(Surv(y1) ~ Surv(y2)))
try(wilcoxsign_test(Surv(y1) ~ y2))
try(wilcoxsign_test(y1 ~ Surv(y2)))
try(friedman_test(Surv(y1) ~ x))
try(quade_test(Surv(y1) ~ x))

### exact two-sample tests with scores gave wrong p-value
y <- 1:20
x <- gl(2, 10)
ox <- ordered(x)

it1 <- independence_test(y ~ x,  distr = exact(algorithm = "shift")) # p = 1e-05
it2 <- independence_test(y ~ x,  distr = exact(algorithm = "shift"), # was p = 1
                         scores = list(x = 1:2))
it3 <- independence_test(y ~ ox, distr = exact(algorithm = "shift")) # was p = 1
it4 <- independence_test(y ~ ox, distr = exact(algorithm = "shift"), # was p = NA
                         scores = list(ox = 3:4))
stopifnot(identical(pvalue(it1), pvalue(it2)))
stopifnot(identical(pvalue(it1), pvalue(it3)))
stopifnot(identical(pvalue(it1), pvalue(it4)))

it5 <- independence_test(y ~ x,  distr = exact(algorithm = "split")) # p = 1e-05
it6 <- independence_test(y ~ x,  distr = exact(algorithm = "split"), # was p = NA
                         scores = list(x = 1:2))
it7 <- independence_test(y ~ ox, distr = exact(algorithm = "split")) # was p = NA
it8 <- independence_test(y ~ ox, distr = exact(algorithm = "split"), # was p = NA
                         scores = list(ox = 3:4))
stopifnot(identical(pvalue(it5), pvalue(it6)))
stopifnot(identical(pvalue(it5), pvalue(it7)))
stopifnot(identical(pvalue(it5), pvalue(it8)))

### 'of_trafo' threw an error for 'x' of length one
of_trafo(gl(3, 1, ordered = TRUE)[1])
