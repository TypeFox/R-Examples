### Regression tests for the r x c x K problem, i.e.,
### testing the independence of a factor
### `y' and a factor factor `x' (possibly blocked)

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)


### generate data: 2 x 2 x K
dat <- data.frame(x = gl(2, 50), y = gl(2, 50)[sample(1:100)],
                  block = gl(10, 10)[sample(1:100)])[sample(1:100, 75),]

### Pearsons Chisq Test, asymptotic distribution
ptwo <- chisq.test(table(dat$x, dat$y), correct = FALSE)$p.value

stopifnot(isequal(pvalue(chisq_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(chisq_test(table(dat$y, dat$x))), ptwo))

### Cochran-Mantel-Haenzel Test, asymptotic distribution
ptwo <- drop(mantelhaen.test(table(dat$x, dat$y, dat$block),
                             correct = FALSE)$p.value)

stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))


### generate data: r x c x K
dat <- data.frame(x = gl(4, 25), y = gl(4, 25)[sample(1:100)],
                  block = gl(2, 50)[sample(1:100)])

### Cochran-Mantel-Haenzel Test, asymptotic distribution
### (was wrong in R < 2.1.0)
ptwo <- drop(mantelhaen.test(table(dat$y, dat$x, dat$block),
                             correct = FALSE)$p.value)

stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))


### generate data: r x c x K
dat <- data.frame(x = gl(4, 25), y = gl(5, 20)[sample(1:100)],
                  block = gl(2, 50)[sample(1:100)])

### Cochran-Mantel-Haenzel Test, asymptotic distribution
### (was wrong in R < 2.1.0)
ptwo <- drop(mantelhaen.test(table(dat$y, dat$x, dat$block),
                             correct = FALSE)$p.value)

stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))


### 2x2 table and maxstat
x <- c(rep(1,51), rep(2,49))
y <- factor(c(rep(0,49), rep(1,51)))[sample(1:100)]
stopifnot(isequal(as.vector(statistic(independence_test(table(x, y)))),
as.vector(statistic(maxstat_test(y ~ x )))))


### maxstat for multiple, ordered and unordered covariates
dat <- data.frame(w = rnorm(100), x = runif(100), y = gl(4, 25)[sample(1:100)],
                  z = ordered(gl(4, 25)[sample(1:100)]))

mt <- maxstat_test(w ~ x, data = dat)
mt
est <- mt@estimates$estimate$cutpoint
stopifnot(isequal(statistic(mt),
                  abs(statistic(independence_test(w ~ (x <= est), data = dat)))))

mt <- maxstat_test(w ~ y, data = dat)
mt
est <- mt@estimates$estimate$cutpoint
xx <- dat$y %in% est
stopifnot(isequal(statistic(mt),
                  abs(statistic(independence_test(w ~ xx, data = dat)))))

mt <- maxstat_test(w ~ z, data = dat)
mt
est <- mt@estimates$estimate$cutpoint
xx <- dat$z <= est
stopifnot(isequal(statistic(mt),
                  abs(statistic(independence_test(w ~ xx, data = dat)))))

mt <- maxstat_test(w ~ x + y + z, data = dat)
mt
est <- mt@estimates$estimate
xsel <- dat[[est[[1]]]]
if (is.factor(xsel) && !is.ordered(xsel)) {
    xx <- xsel %in% est[2]
} else {
    xx <- xsel <= est[2]
}
stopifnot(isequal(statistic(mt),
                  abs(statistic(independence_test(w ~ xx, data = dat)))))


### marginal homogeneity
rating <- c("low", "moderate", "high")
x <- as.table(matrix(c(20, 10,  5,
                       3, 30, 15,
                       0,  5, 40),
                     ncol = 3, byrow = TRUE,
                     dimnames = list(Rater1 = rating, Rater2 = rating)))
### test statistic W_0 = 13.76
### see http://ourworld.compuserve.com/homepages/jsuebersax/mcnemar.htm
stopifnot(all.equal(round(statistic(mh_test(x)), 2), 13.76))
