
### Regression tests for multiple adjustments

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)

### example from Westfall & Wolfinger (1997), Table 4
tab <- as.table(matrix(c(12, 1, 3, 8, 17, 12, 9, 9, 16, 24), nrow = 2,
                       dimnames = list(group = c("Placebo", "Active"),
                                       response = c("Very Poor", "Poor", "Fair",
                                                    "Good", "Excellent"))))
df <- coin:::table2df(tab)

it <- independence_test(response ~ group, data = df,
                        distribution = approximate(B = 100000))

### Table 5, first column: OK
pvalue(it, method = "unadjusted")

### Table 5, next-to-last column: OK
pvalue(it, method = "step-down", distribution = "marginal", type = "Sidak")

### Table 5, last column: OK
pvalue(it, method = "step-down")

### example from Westfall & Wolfinger (1997), Table 1
df <- data.frame(group = factor(c(rep("Control", 50), rep("Treatment", 48))),
                 V1 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 5), rep(1, 5)),
                 V2 = c(rep(0, 50 - 4), rep(1, 4), rep(0, 48 - 3), rep(1, 3)),
                 V3 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 4), rep(1, 4)),
                 V4 = c(rep(0, 50 - 6), rep(1, 6), rep(0, 48 - 4), rep(1, 4)))

### alternative: less
it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df,
                        distribution = approximate(B = 100000),
                        alternative = "less")

### page 4, 2nd column: adjusted p-value = 0.03665 for V1
pvalue(it, method = "single-step", distribution = "marginal", type = "Sidak")

### page 4, 2nd column: adjusted p-value = 0.03698 for V1
### Note: 0.02521 + 0.00532 + 0 + 0.00645 = 0.03698
pvalue(it, method = "single-step", distribution = "marginal")

### alternative: two.sided
it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df,
                        distribution = approximate(B = 100000))

### page 5, 1st column: adjusted p-value = 0.05261 for V1
pvalue(it, method = "single-step", distribution = "marginal", type = "Sidak")

### page 5, 1st column: adjusted p-value = 0.05352 for V1
### Note: 0.02521 + 0.01254 + 0 + 0.01577 = 0.05352
pvalue(it, method = "single-step", distribution = "marginal")

### artificial example, checked against `multtest:mt.maxT'

set.seed(290875)

gr <- gl(2, 50)
x1 <- rnorm(100) + (as.numeric(gr) - 1) * 0.5
x2 <- rnorm(100) - (as.numeric(gr) - 1) * 0.5

pvalue(independence_test(x1 + x2 ~ gr, alternative = "two.sided"),
       method = "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "less"),
       method = "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "greater"),
       method = "single-step")

pvalue(independence_test(x1 + x2 ~ gr, alternative = "two.sided"),
       method = "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "less"),
       method = "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "greater"),
       method = "step-down")

pvalue(independence_test(x1 + x2 ~ gr, alternative = "two.sided",
                         distribution = approximate(B = 10000)),
       method = "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "less",
                         distribution = approximate(B = 10000)),
       method = "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "greater",
                         distribution = approximate(B = 10000)),
       method = "single-step")

pvalue(independence_test(x1 + x2 ~ gr, alternative = "two.sided",
                         distribution = approximate(B = 10000)),
       method = "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "less",
                         distribution = approximate(B = 10000)),
       method = "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alternative = "greater",
                         distribution = approximate(B = 10000)),
       method = "step-down")

if (FALSE) {
    #library("multtest")
    #a <- mt.maxT(t(cbind(x1, x2)), as.numeric(gr) - 1)
    #a[order(a$index),]
    #a <- mt.maxT(t(cbind(x1, x2)), as.numeric(gr) - 1, side = "upper")
    #a[order(a$index),]
    #a <- mt.maxT(t(cbind(x1, x2)), as.numeric(gr) - 1, side = "lower")
    #a[order(a$index),]
}

### Monte Carlo distribution

y <- rnorm(20)
x <- runif(20)

mt <- maxstat_test(y ~ x, distribution = approximate())
pvalue(mt)
pperm(mt, 1)
qperm(mt, 0.9)
dperm(mt, qperm(mt, 0.9))
#support(mt)

mt <- maxstat_test(y ~ x, distribution = approximate(), alternative = "greater")
pvalue(mt)
pperm(mt, 1)
qperm(mt, 0.9)
dperm(mt, qperm(mt, 0.9))
#support(mt)

mt <- maxstat_test(y ~ x, distribution = approximate(), alternative = "less")
pvalue(mt)
pperm(mt, 1)
qperm(mt, 0.9)
dperm(mt, qperm(mt, 0.9))
#support(mt)

### unadjusted

set.seed(290875)

gr <- gl(3, 50)
x1 <- rnorm(150) + (as.numeric(gr) - 1) * 0.5
x2 <- rnorm(150) - (as.numeric(gr) - 1) * 0.5

pvalue(it1 <- independence_test(x1 + x2 ~ gr),
       method = "unadjusted")
pvalue(it2 <- independence_test(x1 + x2 ~ gr, alternative = "less"),
       method = "unadjusted")
pvalue(it3 <- independence_test(x1 + x2 ~ gr, alternative = "greater"),
       method = "unadjusted")

pvalue(it4 <- independence_test(x1 + x2 ~ gr,
                                distribution = approximate(B = 10000)),
       method = "unadjusted")
pvalue(it5 <- independence_test(x1 + x2 ~ gr, alternative = "less",
                                distribution = approximate(B = 10000)),
       method = "unadjusted")
pvalue(it6 <- independence_test(x1 + x2 ~ gr, alternative = "greater",
                                distribution = approximate(B = 10000)),
       method = "unadjusted")

### consistency of minimum p-value for "global"/"single-step"/"step-down"

set.seed(290875); pg1 <- pvalue(it1)[1]
set.seed(290875); pss1 <- pvalue(it1, method = "single-step")
set.seed(290875); psd1 <- pvalue(it1, method = "step-down")
identical(pg1, min(pss1))
identical(pg1, min(psd1))

set.seed(290875); pg2 <- pvalue(it2)[1]
set.seed(290875); pss2 <- pvalue(it2, method = "single-step")
set.seed(290875); psd2 <- pvalue(it2, method = "step-down")
identical(pg2, min(pss2))
identical(pg2, min(psd2))

set.seed(290875); pg3 <- pvalue(it3)[1]
set.seed(290875); pss3 <- pvalue(it3, method = "single-step")
set.seed(290875); psd3 <- pvalue(it3, method = "step-down")
identical(pg3, min(pss3))
identical(pg3, min(psd3))

pg4 <- pvalue(it4)[1]
pss4 <- pvalue(it4, method = "single-step")
psd4 <- pvalue(it4, method = "step-down")
identical(pg4, min(pss4))
identical(pg4, min(psd4))

pg5 <- pvalue(it5)[1]
pss5 <- pvalue(it5, method = "single-step")
psd5 <- pvalue(it5, method = "step-down")
identical(pg5, min(pss5))
identical(pg5, min(psd5))

pg6 <- pvalue(it6)[1]
pss6 <- pvalue(it6, method = "single-step")
psd6 <- pvalue(it6, method = "step-down")
identical(pg6, min(pss6))
identical(pg6, min(psd6))

### adjusted marginal asymptotic p-values

pvalue(it1, method = "single-step", distribution = "marginal")
pvalue(it1, method = "single-step", distribution = "marginal", type = "Sidak")
pvalue(it1, method = "step-down", distribution = "marginal")
pvalue(it1, method = "step-down", distribution = "marginal", type = "Sidak")

pvalue(it2, method = "single-step", distribution = "marginal")
pvalue(it2, method = "single-step", distribution = "marginal", type = "Sidak")
pvalue(it2, method = "step-down", distribution = "marginal")
pvalue(it2, method = "step-down", distribution = "marginal", type = "Sidak")

pvalue(it3, method = "single-step", distribution = "marginal")
pvalue(it3, method = "single-step", distribution = "marginal", type = "Sidak")
pvalue(it3, method = "step-down", distribution = "marginal")
pvalue(it3, method = "step-down", distribution = "marginal", type = "Sidak")

## ### mcp

## YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
##                              42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
##                              38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
##                              31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
##                   site = factor(c(rep("I", 10), rep("II", 10),
##                                   rep("III", 10), rep("IV", 10))))

## ### permutation based Dunnett
## it <- independence_test(length ~ site, data = YOY,
##                         xtrafo = mcp_trafo(site = "Dunnett"),
##                         distribution = approximate(10000),
##                         alternative = "two.sided")
## pvalue(it, method = "npmcp")

## ### asymptotic Dunnett
## it <- independence_test(length ~ site, data = YOY,
##                         xtrafo = mcp_trafo(site = "Dunnett"),
##                         alternative = "two.sided")
## pvalue(it, method = "npmcp")

## ### asymptotic Dunnett, user-defined w/o column names
## cm <- rbind("II  vs I" = c(-1, 1, 0, 0),
##             "III vs I" = c(-1, 0, 1, 0),
##             "IV  vs I" = c(-1, 0, 0, 1))
## it <- independence_test(length ~ site, data = YOY,
##                         xtrafo = mcp_trafo(site = cm),
##                         alternative = "two.sided")
## pvalue(it, method = "npmcp")
