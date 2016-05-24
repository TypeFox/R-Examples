### Regression tests for the distribution functions

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)


### generate independent two-sample data
dta <- data.frame(y = rnorm(20), x = gl(2, 10), b = factor(rep(1:4, 5)),
                  w = rep(1:3, length = 20))
dta$y5 <- round(dta$y, 5)
dta$y3 <- round(dta$y, 3)

### generate paired two-sample data
dta2 <- data.frame(y = as.vector(rbind(abs(dta$y) * (dta$y >= 0),
                                       abs(dta$y) * (dta$y < 0))),
                   x = factor(rep(0:1, length(dta$y)),
                              labels = c("pos", "neg")),
                   b = gl(length(dta$y), 2))
dta2$y5 <- round(dta2$y, 5)
dta2$y3 <- round(dta2$y, 3)


### check 'algorithm = "auto"'

### scores that can't be mapped into integers

### two-sample with block
try(independence_test(y ~ x | b, data = dta,
                      distribution = exact(algo = "auto")))

try(independence_test(y ~ x | b, data = dta,
                      distribution = exact(algo = "shift")))

try(independence_test(y ~ x | b, data = dta,
                      distribution = exact(algo = "split-up")))

### two-sample without block
it <- independence_test(y ~ x, data = dta,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

try(independence_test(y ~ x, data = dta,
                      distribution = exact(algo = "shift")))

it <- independence_test(y ~ x, data = dta,
                        distribution = exact(algo = "split-up"))
it@distribution@name
pvalue(it)

### paired two-sample
try(symmetry_test(y ~ x | b, data = dta2, paired = TRUE,
                  distribution = exact(algo = "auto")))

try(symmetry_test(y ~ x | b, data = dta2, paired = TRUE,
                  distribution = exact(algo = "shift")))

try(symmetry_test(y ~ x | b, data = dta2, paired = TRUE,
                  distribution = exact(algo = "split-up")))

### mapped into integers using 'fact'

### two-sample with block
it <- independence_test(y5 ~ x | b, data = dta,
                        distribution = exact(algo = "auto", fact = 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x | b, data = dta,
                        distribution = exact(algo = "shift", fact = 1e5))
it@distribution@name
pvalue(it)

try(independence_test(y5 ~ x | b, data = dta,
                        distribution = exact(algo = "split-up")))

### two-sample without block
it <- independence_test(y5 ~ x, data = dta,
                        distribution = exact(algo = "auto", fact = 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x, data = dta,
                        distribution = exact(algo = "shift", fact = 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x, data = dta,
                        distribution = exact(algo = "split-up"))
it@distribution@name
pvalue(it)

### paired two-sample
st <- symmetry_test(y5 ~ x | b, data = dta2, paired = TRUE,
                    distribution = exact(algo = "auto", fact = 1e5))
st@distribution@name
pvalue(st)

st <- symmetry_test(y5 ~ x | b, data = dta2, paired = TRUE,
                    distribution = exact(algo = "shift", fact = 1e5))
st@distribution@name
pvalue(st)

try(symmetry_test(y5 ~ x | b, data = dta2, paired = TRUE,
                  distribution = exact(algo = "split-up")))

### automatically mapped into integers

### two-sample with block
it <- independence_test(y3 ~ x | b, data = dta,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x | b, data = dta,
                        distribution = exact(algo = "shift"))
it@distribution@name
pvalue(it)

try(independence_test(y3 ~ x | b, data = dta,
                      distribution = exact(algo = "split-up")))

### two-sample without block
it <- independence_test(y3 ~ x, data = dta,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x, data = dta,
                        distribution = exact(algo = "shift"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x, data = dta,
                        distribution = exact(algo = "split-up"))
it@distribution@name
pvalue(it)

### paired two-sample
st <- symmetry_test(y3 ~ x | b, data = dta2, paired = TRUE,
                    distribution = exact(algo = "auto"))
st@distribution@name
pvalue(st)

st <- symmetry_test(y3 ~ x | b, data = dta2, paired = TRUE,
                    distribution = exact(algo = "shift"))
st@distribution@name
pvalue(st)

try(symmetry_test(y3 ~ x | b, data = dta2, paired = TRUE,
                  distribution = exact(algo = "split-up")))


### check exact tests with weights
itw1 <- independence_test(y3 ~ x, data = dta, weights = ~ w,
                          distribution = exact(algorithm = "shift"))
itw2 <- independence_test(y3 ~ x, data = dta, weights = ~ w,
                          distribution = exact(algorithm = "split-up"))
y3w <- with(dta, rep(y3, w))
xw <- with(dta, rep(x, w))
it1 <- independence_test(y3w ~ xw, distribution = exact(algorithm = "shift"))
it2 <- independence_test(y3w ~ xw, distribution = exact(algorithm = "split-up"))
stopifnot(isequal(pvalue(itw1), pvalue(it1)))
stopifnot(isequal(pvalue(itw1), pvalue(it2)))
stopifnot(isequal(pvalue(itw2), pvalue(it1)))
stopifnot(isequal(pvalue(itw2), pvalue(it2)))

Convictions <-
    matrix(c(2, 10, 15, 3), nrow = 2,
           dimnames = list(c("Dizygotic", "Monozygotic"),
                           c("Convicted", "Not convicted")))
itw1 <- independence_test(as.table(Convictions), alternative = "less",
                          distribution = exact(algorithm = "shift"))
itw2 <- independence_test(as.table(Convictions), alternative = "less",
                          distribution = exact(algorithm = "split-up"))
it1 <- independence_test(Var2 ~ Var1, alternative = "less",
                         data = coin:::table2df(as.table(Convictions)),
                         distribution = exact(algorithm = "shift"))
it2 <- independence_test(Var2 ~ Var1, alternative = "less",
                         data = coin:::table2df(as.table(Convictions)),
                         distribution = exact(algorithm = "split-up"))
stopifnot(isequal(pvalue(itw1), pvalue(it1)))
stopifnot(isequal(pvalue(itw1), pvalue(it2)))
stopifnot(isequal(pvalue(itw2), pvalue(it1)))
stopifnot(isequal(pvalue(itw2), pvalue(it2)))


### check support, pperm, dperm, qperm,
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
dta1 <- data.frame(
    y = c(y1, y2),
    x = gl(2, length(y1)),
    b = factor(rep(seq_along(y1), 2)))

diff <- y1 - y2
dta2 <- data.frame(
    y = as.vector(rbind(abs(diff) * (diff >= 0), abs(diff) * (diff < 0))),
    x = factor(rep(0:1, length(diff)), labels = c("pos", "neg")),
    b <- gl(length(diff), 2))

### shift without block
it1_SR <- independence_test(y ~ x, data = dta1,
                            distribution = exact(algorithm = "shift"))
supp_it1_SR <- support(it1_SR)
stopifnot(!is.unsorted(supp_it1_SR))
stopifnot(all(supp_it1_SR == unique(supp_it1_SR)))
round(pp_it1_SR <- pperm(it1_SR, supp_it1_SR), 7)
round(dp_it1_SR <- dperm(it1_SR, supp_it1_SR), 7)
round(qp_it1_SR <- qperm(it1_SR, seq(0, 1, 0.01)), 7)

### split-up without block
it1_vdW <- independence_test(y ~ x, data = dta1,
                             distribution = exact(algorithm = "split-up"))
round(pp_it1_vdW <- pperm(it1_vdW, supp_it1_SR), 7)
round(qp_it1_vdW <- qperm(it1_vdW, seq(0, 1, 0.01)), 7)

### should be equal
stopifnot(isequal(pp_it1_SR, pp_it1_vdW))
## <FIXME> Doesn't pass under Solaris or Linux w/o long doubles
## stopifnot(isequal(qp_it1_SR[-c(1, 101)], qp_it1_vdW[-c(1, 101)]))
## </FIXME>
stopifnot(isequal(pvalue(it1_SR), pvalue(it1_vdW)))

### shift with block
it2_SR <- independence_test(y ~ x | b, data = dta1,
                            distribution = exact(algorithm = "shift"))
supp_it2_SR <- support(it2_SR)
stopifnot(!is.unsorted(supp_it2_SR))                         # failed in < 1.1-0
stopifnot(all(supp_it2_SR == unique(supp_it2_SR)))           # failed in < 1.1-0

round(pp_it2_SR <- pperm(it2_SR, supp_it2_SR), 7)
round(dp_it2_SR <- dperm(it2_SR, supp_it2_SR), 7)            # failed in < 1.1-0
round(qp_it2_SR <- qperm(it2_SR, seq(0, 1, 0.01)), 7)

### paired shift with block
st3_SR <- symmetry_test(y ~ x | b, data = dta2,
                        distribution = exact(algorithm = "shift"),
                        paired = TRUE)
supp_st3_SR <- support(st3_SR)
stopifnot(!is.unsorted(supp_st3_SR))
stopifnot(all(supp_st3_SR == unique(supp_st3_SR)))
round(pp_st3_SR <- pperm(st3_SR, supp_st3_SR), 7)
round(dp_st3_SR <- dperm(st3_SR, supp_st3_SR), 7)
round(qp_st3_SR <- qperm(st3_SR, seq(0, 1, 0.01)), 7)

### should be equal
stopifnot(isequal(pp_it2_SR, pp_st3_SR))                     # failed in < 1.1-0
stopifnot(isequal(qp_it2_SR, qp_st3_SR))                     # failed in < 1.1-0
stopifnot(isequal(pvalue(it2_SR), pvalue(st3_SR)))


### exact test based on quadratic forms

### shift with block
itq1 <- independence_test(y ~ x | b, data = dta1,
                          distribution = exact(algorithm = "shift"),
                          teststat = "quad")
its1 <- independence_test(y ~ x | b, data = dta1,
                          distribution = exact(algorithm = "shift"),
                          teststat = "scalar")
stopifnot(isequal(statistic(itq1), statistic(its1)^2))
stopifnot(isequal(pvalue(itq1), pvalue(its1)))
stopifnot(isequal(support(itq1), support(its1)[support(its1) >= 0]^2))

### paired shift with block
sts2 <- symmetry_test(y ~ x | b, data = dta2,
                      distribution = exact(algorithm = "shift"),
                      paired = TRUE, teststat = "scalar")
stq2 <- symmetry_test(y ~ x | b, data = dta2,
                      distribution = exact(algorithm = "shift"),
                      paired = TRUE, teststat = "quad")
stopifnot(isequal(statistic(stq2), statistic(sts2)^2))
stopifnot(isequal(pvalue(stq2), pvalue(sts2)))
stopifnot(isequal(support(stq2), support(sts2)[support(sts2) >= 0]^2))

### should be equal
stopifnot(isequal(statistic(itq1), statistic(stq2)))
stopifnot(isequal(pvalue(itq1), pvalue(stq2)))
stopifnot(isequal(support(itq1), support(stq2)))

### dperm gave an error here (fixed in r869)
tab <- as.table(matrix(c(3, 0, 0, 3), ncol = 2))
itq <- independence_test(tab, distribution = exact(algorithm = "shift"),
                         teststat = "quad")
stopifnot(is.numeric(dperm(itq, statistic(itq))))


### check vectorization
dta <- data.frame(y1 = sample(1:20), y2 = sample(1:20), x = gl(2, 10))

### univariate, asymptotic and scalar
it_uas <- independence_test(y1 ~ x,  data = dta,
                            distribution = "asymptotic", teststat = "scalar")
round(dp_uas <- dperm(it_uas, 0:1), 7)
stopifnot(isequal(dp_uas, c(dperm(it_uas, 0), dperm(it_uas, 1))))
round(pp_uas <- pperm(it_uas, 0:1), 7)
stopifnot(isequal(pp_uas, c(pperm(it_uas, 0), pperm(it_uas, 1))))
round(qp_uas <- qperm(it_uas, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_uas, c(qperm(it_uas, 0.5), qperm(it_uas, 0.75))))

### univariate, asymptotic and quad
it_uaq <- independence_test(y1 ~ x,  data = dta,
                            distribution = "asymptotic", teststat = "quad")
round(dp_uaq <- dperm(it_uaq, 0:1), 7)
stopifnot(isequal(dp_uaq, c(dperm(it_uaq, 0), dperm(it_uaq, 1))))
round(pp_uaq <- pperm(it_uaq, 0:1), 7)
stopifnot(isequal(pp_uaq, c(pperm(it_uaq, 0), pperm(it_uaq, 1))))
round(qp_uaq <- qperm(it_uaq, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_uaq, c(qperm(it_uaq, 0.5), qperm(it_uaq, 0.75))))

### multivariate, asymptotic and max
it_mam <- independence_test(y1 + y2 ~ x,  data = dta,
                            distribution = "asymptotic", teststat = "max")
round(dp_mam <- dperm(it_mam, 0:1), 7)
#stopifnot(isequal(dp_mam, c(dperm(it_mam, 1), dperm(it_mam, 2))))
round(pp_mam <- pperm(it_mam, 0:1), 7)                       # failed in < 1.1-0
stopifnot(isequal(pp_mam, c(pperm(it_mam, 0), pperm(it_mam, 1))))
round(qp_mam <- qperm(it_mam, c(0.5, 0.75)), 7)              # failed in < 1.1-0
stopifnot(isequal(qp_mam, c(qperm(it_mam, 0.5), qperm(it_mam, 0.75))))

### multivariate, asymptotic and quad
it_maq <- independence_test(y1 + y2 ~ x,  data = dta,
                            distribution = "asymptotic", teststat = "quad")
round(dp_maq <- dperm(it_maq, 0:1), 7)
stopifnot(isequal(dp_maq, c(dperm(it_maq, 0), dperm(it_maq, 1))))
round(pp_maq <- pperm(it_maq, 0:1), 7)
stopifnot(isequal(pp_maq, c(pperm(it_maq, 0), pperm(it_maq, 1))))
round(qp_maq <- qperm(it_maq, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_maq, c(qperm(it_maq, 0.5), qperm(it_maq, 0.75))))

### univariate, approximate and scalar
it_ums <- independence_test(y1 ~ x,  data = dta,
                            distribution = "approximate", teststat = "scalar")
round(dp_ums <- dperm(it_ums, 0:1), 7)
stopifnot(isequal(dp_ums, c(dperm(it_ums, 0), dperm(it_ums, 1))))
round(pp_ums <- pperm(it_ums, 0:1), 7)
stopifnot(isequal(pp_ums, c(pperm(it_ums, 0), pperm(it_ums, 1))))
round(qp_ums <- qperm(it_ums, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_ums, c(qperm(it_ums, 0.5), qperm(it_ums, 0.75))))

### univariate, approximate and quad
it_umq <- independence_test(y1 ~ x,  data = dta,
                            distribution = "approximate", teststat = "quad")
round(dp_umq <- dperm(it_umq, 0:1), 7)
stopifnot(isequal(dp_umq, c(dperm(it_umq, 0), dperm(it_umq, 1))))
round(pp_umq <- pperm(it_umq, 0:1), 7)
stopifnot(isequal(pp_umq, c(pperm(it_umq, 0), pperm(it_umq, 1))))
round(qp_umq <- qperm(it_umq, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_umq, c(qperm(it_umq, 0.5), qperm(it_umq, 0.75))))

### multivariate, approximate and max
it_mmm <- independence_test(y1 + y2 ~ x,  data = dta,
                            distribution = "approximate", teststat = "max")
round(dp_mmm <- dperm(it_mmm, c(0, 1)), 7)
stopifnot(isequal(dp_mmm, c(dperm(it_mmm, 0), dperm(it_mmm, 1))))
round(pp_mmm <- pperm(it_mmm, 0:1), 7)
stopifnot(isequal(pp_mmm, c(pperm(it_mmm, 0), pperm(it_mmm, 1))))
round(qp_mmm <- qperm(it_mmm, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_mmm, c(qperm(it_mmm, 0.5), qperm(it_mmm, 0.75))))

### multivariate, approximate and quad
it_mmq <- independence_test(y1 + y2 ~ x,  data = dta,
                            distribution = "approximate", teststat = "quad")
round(dp_mmq <- dperm(it_mmq, 0:1), 7)
stopifnot(isequal(dp_mmq, c(dperm(it_mmq, 0), dperm(it_mmq, 1))))
round(pp_mmq <- pperm(it_mmq, 0:1), 7)
stopifnot(isequal(pp_mmq, c(pperm(it_mmq, 0), pperm(it_mmq, 1))))
round(qp_mmq <- qperm(it_mmq, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_mmq, c(qperm(it_mmq, 0.5), qperm(it_mmq, 0.75))))

### univariate, exact and scalar
it_ues <- independence_test(y1 ~ x,  data = dta,
                            distribution = "exact", teststat = "scalar")
round(dp_ues <- dperm(it_ues, 0:1), 7)
stopifnot(isequal(dp_ues, c(dperm(it_ues, 0), dperm(it_ues, 1))))
round(pp_ues <- pperm(it_ues, 0:1), 7)
stopifnot(isequal(pp_ues, c(pperm(it_ues, 0), pperm(it_ues, 1))))
round(qp_ues <- qperm(it_ues, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_ues, c(qperm(it_ues, 0.5), qperm(it_ues, 0.75))))

### univariate, exact and quad
it_ueq <- independence_test(y1 ~ x,  data = dta,
                            distribution = "exact", teststat = "quad")
round(dp_ueq <- dperm(it_ueq, 0:1), 7)
stopifnot(isequal(dp_ueq, c(dperm(it_ueq, 0), dperm(it_ueq, 1))))
round(pp_ueq <- pperm(it_ueq, 0:1), 7)
stopifnot(isequal(pp_ueq, c(pperm(it_ueq, 0), pperm(it_ueq, 1))))
round(qp_ueq <- qperm(it_ueq, c(0.5, 0.75)), 7)
stopifnot(isequal(qp_ueq, c(qperm(it_ueq, 0.5), qperm(it_ueq, 0.75))))
