
### Regression tests for the 2 sample problem, i.e.,
### testing the independence of a numeric variable
### `y' and a binary factor `x' (possibly blocked)

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)

### generate data
dat <- data.frame(x = gl(2, 50), y = rnorm(100), block = gl(5, 20))[sample(1:100, 75), ]


### Wilcoxon Mann-Whitney Rank Sum Test

### asymptotic distribution
ptwo <- wilcox.test(y ~ x, data = dat, correct = FALSE, exact = FALSE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less",
                     correct = FALSE, exact = FALSE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater",
                        correct = FALSE, exact = FALSE)$p.value

stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "less")),
                  pless))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater")),
                  pgreater))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo))), ptwo))
### check direct supply of a function via ytrafo
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = rank_trafo)), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
    alternative = "less")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
    alternative = "greater")), pgreater))

### exact distribution
ptwo <- wilcox.test(y ~ x, data = dat, exact = TRUE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less", exact = TRUE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater",
                        exact = TRUE)$p.value

stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, distribution = "exact")),
                  ptwo))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "less",
                             distribution = "exact")), pless))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater",
                             distribution = "exact")), pgreater))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo))), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
    alternative = "less")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
    alternative = "greater")), pgreater))

### approximated distribution
rtwo <- pvalue(wilcox_test(y ~ x, data = dat, distribution = "approx")) / ptwo
rless <- pvalue(wilcox_test(y ~ x, data = dat, alternative = "less",
                   distribution = "approx")) / pless
rgreater <- pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater",
                   distribution = "approx")) / pgreater
stopifnot(all(c(rtwo, rless, rgreater) > 0.9 &
              c(rtwo, rless, rgreater) < 1.1))

### <FIXME> add block examples </FIXME>

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact"))

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt",
                   alternative = "less"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx",
                   alternative = "less"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact",
                   alternative = "less"))

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt",
                   alternative = "greater"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx",
                   alternative = "greater"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact",
                   alternative = "greater"))

### sanity checks
try(wilcox_test(x ~ y, data = dat))
try(wilcox_test(x ~ y | y, data = dat))


### Ansari-Bradley Test

### asymptotic distribution
ptwo <- ansari.test(y ~ x, data = dat, correct = FALSE, exact = FALSE)$p.value
pless <- ansari.test(y ~ x, data = dat, alternative = "less",
                     correct = FALSE, exact = FALSE)$p.value
pgreater <- ansari.test(y ~ x, data = dat, alternative = "greater",
                        correct = FALSE, exact = FALSE)$p.value

stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "less")),
                  pless))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "greater")),
                  pgreater))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo))), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo),
    alternative = "greater")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo),
    alternative = "less")), pgreater))

### exact distribution
ptwo <- ansari.test(y ~ x, data = dat, exact = TRUE)$p.value
pless <- ansari.test(y ~ x, data = dat, alternative = "less", exact = TRUE)$p.value
pgreater <- ansari.test(y ~ x, data = dat, alternative = "greater",
                        exact = TRUE)$p.value

### <FIXME>: Definition of two-sided P-values! </FIXME>
(isequal(pvalue(ansari_test(y ~ x, data = dat, distribution = "exact")),
                  ptwo))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "less",
                             distribution = "exact")), pless))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "greater",
                             distribution = "exact")), pgreater))

### <FIXME>: Definition of two-sided P-values! </FIXME>
(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact",
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo))), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact",
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo),
    alternative = "greater")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact",
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo),
    alternative = "less")), pgreater))

### approximated distribution
rtwo <- pvalue(ansari_test(y ~ x, data = dat, distribution = "approx")) / ptwo
rless <- pvalue(ansari_test(y ~ x, data = dat, alternative = "less",
                   distribution = "approx")) / pless
rgreater <- pvalue(ansari_test(y ~ x, data = dat, alternative = "greater",
                   distribution = "approx")) / pgreater
### <FIXME> ??? </FIXME>
(all(c(rtwo, rless, rgreater) > 0.9 &
              c(rtwo, rless, rgreater) < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(ansari_test(x ~ y, data = dat))
try(ansari_test(x ~ y | y, data = dat))


### One-way Test
oneway_test(y ~ x, dat = dat)
oneway_test(y ~ x, dat = dat, alternative = "less")
oneway_test(y ~ x, dat = dat, alternative = "greater")


### Normal Scores Test
normal_test(y ~ x, dat = dat)
normal_test(y ~ x, dat = dat, alternative = "less")
normal_test(y ~ x, dat = dat, alternative = "greater")


### Median Test
median_test(y ~ x, dat = dat)
median_test(y ~ x, dat = dat, alternative = "less")
median_test(y ~ x, dat = dat, alternative = "greater")


### Savage Test
savage_test(y ~ x, dat = dat)
savage_test(y ~ x, dat = dat, alternative = "less")
savage_test(y ~ x, dat = dat, alternative = "greater")


### Taha Test
taha_test(y ~ x, dat = dat)
taha_test(y ~ x, dat = dat, alternative = "less")
taha_test(y ~ x, dat = dat, alternative = "greater")


### Klotz Test
klotz_test(y ~ x, dat = dat)
klotz_test(y ~ x, dat = dat, alternative = "less")
klotz_test(y ~ x, dat = dat, alternative = "greater")


### Mood Test
mood_test(y ~ x, dat = dat)
mood_test(y ~ x, dat = dat, alternative = "less")
mood_test(y ~ x, dat = dat, alternative = "greater")


### Fligner-Killeen Test
fligner_test(y ~ x, dat = dat)
fligner_test(y ~ x, dat = dat, alternative = "less")
fligner_test(y ~ x, dat = dat, alternative = "greater")


### Conover-Iman Test
conover_test(y ~ x, dat = dat)
conover_test(y ~ x, dat = dat, alternative = "less")
conover_test(y ~ x, dat = dat, alternative = "greater")


### Logrank Test
logrank_test(Surv(y) ~ x, dat = dat)
logrank_test(Surv(y) ~ x, dat = dat, alternative = "less")
logrank_test(Surv(y) ~ x, dat = dat, alternative = "greater")


### confidence intervals, cf Bauer 1972

### Location Tests
location <- data.frame(y = c(6, 20, 27, 38, 46, 51, 54, 57,
                             10, 12, 15, 21, 32, 40, 41, 45),
                       x = gl(2, 8))

### Wilcoxon Rank-Sum Test
wt <- wilcox_test(y ~ x, data = location, conf.int = TRUE, di = "ex")
wt
ci <- confint(wt)
wt0 <- wilcox.test(y ~ x, data = location, conf.int = TRUE)
stopifnot(isequal(wt0$confint, ci$confint))
stopifnot(isequal(wt0$estimate, ci$estimate))

### Normal Scores Test
nt <- normal_test(y ~ x, data = location, conf.int = TRUE, di = "ex")
nt
ci <- confint(nt)
stopifnot(isequal(ci$conf.int, c(-6, 30)))
stopifnot(isequal(ci$estimate, 11))

### Median Test
median_test(y ~ x, data = location, conf.int = TRUE, di = "ex")

### Savage Test
savage_test(y ~ x, data = location, conf.int = TRUE, di = "ex")


### Scale Tests
scale <- data.frame(y = c(-101, -35, -13, 10, 130, 236, 370, 556,
                          -145, -140, -40, -30, 2, 27, 68, 290),
                    x = gl(2, 8))

### Ansari-Bradley Test
at <- ansari_test(y ~ x, data = scale, di = "ex",
                  conf.int = TRUE, conf.level = 0.988)
at
ci <- confint(at)
stopifnot(isequal(ci$conf.int, c(10, 556) / c(68, 27)))
stopifnot(isequal(ci$estimate, mean(c(35/30, 370 / 290))))

### Taha Test
taha_test(y ~ x, data = scale, di = "ex",
          conf.int = TRUE, conf.level = 0.54)

### Klotz Test
klotz_test(y ~ x, data = scale, di = "ex",
           conf.int = TRUE, conf.level = 0.988)

### Mood Test
mood_test(y ~ x, data = scale, di = "ex",
          conf.int = TRUE, conf.level = 0.988)

### Fligner-Killeen Test
fligner_test(y ~ x, data = scale, di = "ex",
                  conf.int = TRUE, conf.level = 0.988)

### Conover-Iman Test
conover_test(y ~ x, data = scale, di = "ex",
             conf.int = TRUE, conf.level = 0.988)


### ties handling
y1 <- c(14 , 18 , 2 , 4 , -5 , 14 , -3 , -1 , 1 , 6 , 3 , 3)
x1 <- c(8 , 26 , -7 , -1 , 2 , 9 , 0 , -4 , 13 , 3 , 3 , 4)
pvalue(wilcoxsign_test(y1~x1,alter="greater",dist=exact(),
                       zero.method = "Wilcoxon"))
pvalue(wilcoxsign_test(y1~x1,alter="greater",dist=exact()))


### Weighted logrank tests

### Collett (2003, p. 9, Table 1.3)
prostatic <- data.frame(
    time = c(13, 52,  6, 40, 10,  7, 66, 10, 10, 14,
             16 , 4, 65,  5, 11, 10, 15,  5, 76, 56,
             88, 24, 51,  4, 40,  8, 18,  5, 16, 50,
             40,  1, 36,  5, 10, 91, 18,  1, 18,  6,
              1, 23, 15, 18, 12, 12, 17,  3),
    event = c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1,
              1, 1, 1, 1, 0, 1, 0, 1, 0, 0,
              1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
              1, 1, 1, 1, 0, 1, 1, 0),
    Hb = c(14.6, 12.0, 11.4, 10.2, 13.2,  9.9, 12.8, 14.0,  7.5, 10.6,
           11.2, 10.1,  6.6,  9.7,  8.8,  9.6, 13.0, 10.4, 14.0, 12.5,
           14.0, 12.4, 10.1,  6.5, 12.8,  8.2, 14.4, 10.2, 10.0,  7.7,
            5.0,  9.4, 11.0,  9.0, 14.0, 11.0, 10.8,  5.1, 13.0,  5.1,
           11.3, 14.6,  8.8,  7.5,  4.9,  5.5,  7.5, 10.2))
prostatic <- within(prostatic,
                    group <- factor(Hb > 11.0, labels = as.roman(1:2)))

### Leton and Zuluaga (2005, p. 384, Table 9)

### Gehan
lt <- logrank_test(Surv(time, event) ~ group, data = prostatic,
                   type = "Gehan")
isequal(round(statistic(lt)^2, 4), 3.8400)
isequal(round(pvalue(lt), 4), 0.0500)

### Peto-Peto
lt <- logrank_test(Surv(time, event) ~ group, data = prostatic,
                   type = "Fleming-Harrington", rho = 1)
isequal(round(statistic(lt)^2, 4), 4.0657)
isequal(round(pvalue(lt), 4), 0.0438)

### Prentice
lt <- logrank_test(Surv(time, event) ~ group, data = prostatic,
                   type = "Prentice")
isequal(round(statistic(lt)^2, 4), 4.1229)
isequal(round(pvalue(lt), 4), 0.0423)

### LR Altshuler
lt <- logrank_test(Surv(time, event) ~ group, data = prostatic)
isequal(round(statistic(lt)^2, 4), 4.4343)
isequal(round(pvalue(lt), 4), 0.0352)

### Tarone-Ware
lt <- logrank_test(Surv(time, event) ~ group, data = prostatic,
                   type = "Tarone-Ware")
isequal(round(statistic(lt)^2, 4), 4.3443)
isequal(round(pvalue(lt), 4), 0.0371)


### Paired tests

### sanity check
set.seed(123)
x <- factor(rep(1:2, 15))
y <- as.integer(round((rnorm(30) + as.numeric(x)) * 1000))
id <- gl(15, 2)
try(symmetry_test(y ~ x | id, distribution = "exact", paired = TRUE))
