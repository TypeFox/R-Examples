
### Regression tests for the K sample problem, i.e.,
### testing the independence of a numeric variable
### `y' and a factor `x' (possibly blocked)

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)

### generate data
dat <- data.frame(x = gl(4, 25), y = rnorm(100), block = gl(5, 20))[sample(1:100, 50), ]


### Kruskal-Wallis Test

### asymptotic distribution
ptwo <- kruskal.test(y ~ x, data = dat)$p.value

stopifnot(isequal(pvalue(kruskal_test(y ~ x, data = dat)), ptwo))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo))), ptwo))

### approximated distribution
rtwo <- pvalue(kruskal_test(y ~ x, data = dat, distribution = "approx")) / ptwo
stopifnot(all(rtwo > 0.9 &
              rtwo < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(kruskal_test(x ~ y, data = dat))
try(kruskal_test(x ~ y | y, data = dat))


### Fligner-Killeen Test

### asymptotic distribution
ptwo <- fligner.test(y ~ x, data = dat)$p.value

stopifnot(isequal(pvalue(fligner_test(y ~ x, data = dat)), ptwo))

dat$yy <- dat$y - tapply(dat$y, dat$x, median)[dat$x]
stopifnot(isequal(pvalue(oneway_test(yy ~ x, data = dat, distribution = "asympt",
    ytrafo = function(data) trafo(data, numeric_trafo = fligner_trafo))), ptwo))

### approximated distribution
rtwo <- pvalue(fligner_test(y ~ x, data = dat, distribution = "approx")) / ptwo
stopifnot(all(rtwo > 0.9 &
              rtwo < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(fligner_test(x ~ y, data = dat))
try(fligner_test(x ~ y | y, data = dat))


### One-way Test
oneway_test(y ~ x, data = dat)

oneway_test(y ~ ordered(x), data = dat)
oneway_test(y ~ ordered(x), data = dat,
            alternative = "less")
oneway_test(y ~ ordered(x), data = dat,
            alternative = "greater")

oneway_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
oneway_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
oneway_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Normal Scores Test
normal_test(y ~ x, data = dat)

normal_test(y ~ ordered(x), data = dat)
normal_test(y ~ ordered(x), data = dat,
            alternative = "less")
normal_test(y ~ ordered(x), data = dat,
            alternative = "greater")

normal_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
normal_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
normal_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Median Test
median_test(y ~ x, data = dat)

median_test(y ~ ordered(x), data = dat)
median_test(y ~ ordered(x), data = dat,
            alternative = "less")
median_test(y ~ ordered(x), data = dat,
            alternative = "greater")

median_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
median_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
median_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Savage Test
savage_test(y ~ x, data = dat)

savage_test(y ~ ordered(x), data = dat)
savage_test(y ~ ordered(x), data = dat,
            alternative = "less")
savage_test(y ~ ordered(x), data = dat,
            alternative = "greater")

savage_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
savage_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
savage_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Taha Test
taha_test(y ~ x, data = dat)

try(taha_test(y ~ ordered(x), data = dat))
try(taha_test(y ~ ordered(x), data = dat,
              alternative = "less"))
try(taha_test(y ~ ordered(x), data = dat,
              alternative = "greater"))

try(taha_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8))))
try(taha_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
              alternative = "less"))
try(taha_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
              alternative = "greater"))


### Klotz Test
klotz_test(y ~ x, data = dat)

try(klotz_test(y ~ ordered(x), data = dat))
try(klotz_test(y ~ ordered(x), data = dat,
               alternative = "less"))
try(klotz_test(y ~ ordered(x), data = dat,
               alternative = "greater"))

try(klotz_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8))))
try(klotz_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
               alternative = "less"))
try(klotz_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
               alternative = "greater"))


### Mood Test
mood_test(y ~ x, data = dat)

try(mood_test(y ~ ordered(x), data = dat))
try(mood_test(y ~ ordered(x), data = dat,
              alternative = "less"))
try(mood_test(y ~ ordered(x), data = dat,
              alternative = "greater"))

try(mood_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8))))
try(mood_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
              alternative = "less"))
try(mood_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
              alternative = "greater"))


### Ansari-Bradley Test
ansari_test(y ~ x, data = dat)

try(ansari_test(y ~ ordered(x), data = dat))
try(ansari_test(y ~ ordered(x), data = dat,
                alternative = "less"))
try(ansari_test(y ~ ordered(x), data = dat,
                alternative = "greater"))

try(ansari_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8))))
try(ansari_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
                alternative = "less"))
try(ansari_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
                alternative = "greater"))


### Conover-Iman Test
conover_test(y ~ x, data = dat)

try(conover_test(y ~ ordered(x), data = dat))
try(conover_test(y ~ ordered(x), data = dat,
                 alternative = "less"))
try(conover_test(y ~ ordered(x), data = dat,
                 alternative = "greater"))

try(conover_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8))))
try(conover_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
                 alternative = "less"))
try(conover_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
                 alternative = "greater"))


### Logrank Test
logrank_test(Surv(y) ~ x, data = dat)

logrank_test(Surv(y) ~ ordered(x), data = dat)
logrank_test(Surv(y) ~ ordered(x), data = dat,
             alternative = "less")
logrank_test(Surv(y) ~ ordered(x), data = dat,
             alternative = "greater")

logrank_test(Surv(y) ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
logrank_test(Surv(y) ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
             alternative = "less")
logrank_test(Surv(y) ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
             alternative = "greater")


### Weighted logrank tests

### Lee & Wang (2003, p. 130, Table 5.11)
leukemia <- data.frame(
    time = c( 4,   5,   9,  10,  12,  13, 10,
             23,  28,  28,  28,  29,
             31,  32,  37,  41,  41,
             57,  62,  74, 100, 139,
             20, 258, 269,
              8,  10,  10,  12,  14,
             20,  48,  70,  75,  99, 103,
            162, 169, 195, 220,
            161, 199, 217,
            245,
              8,  10,  11,  23,  25,  25,
             28,  28,  31,  31,  40,
             48,  89, 124, 143,
             12, 159, 190, 196,
            197, 205, 219),
    event = c(1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1,
              1, 1, 1, 1, 1,
              1, 1, 1, 1, 1,
              0, 0, 0,
              1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1,
              0, 0, 0,
              0,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1,
              1, 1, 1, 1,
              0, 0, 0, 0,
              0, 0, 0),
    group = factor(rep(1:3, c(25, 19, 22)), labels = as.roman(1:3)))

### Leton and Zuluaga (2002, p. 25, Table 6)

### Gehan, X^2_SG
lt <- logrank_test(Surv(time, event) ~ group, data = leukemia,
                   type = "Gehan")
stopifnot(all(-statistic(lt, "linear") == c(273, -170, -103)))
isequal(round(statistic(lt), 3), 3.612)
isequal(round(pvalue(lt), 4), 0.1643)

### Peto-Peto, X^2_SPP
lt <- logrank_test(Surv(time, event) ~ group, data = leukemia,
                   type = "Fleming-Harrington", rho = 1)
stopifnot(all(round(-statistic(lt, "linear"), 3) == c(4.171, -2.582, -1.589)))
isequal(round(statistic(lt), 3), 3.527)
isequal(round(pvalue(lt), 4), 0.1715)

### X^2_S1
lt <- logrank_test(Surv(time, event) ~ group, data = leukemia,
                   type = "Prentice")
stopifnot(all(round(-statistic(lt, "linear"), 3) == c(4.100, -2.503, -1.597)))
isequal(round(statistic(lt), 3), 3.639)
isequal(round(pvalue(lt), 4), 0.1621)

### LR Altshuler, X^2_SLRA
lt <- logrank_test(Surv(time, event) ~ group, data = leukemia)
stopifnot(all(round(-statistic(lt, "linear"), 3) == c(6.635, -3.693, -2.942)))
isequal(round(statistic(lt), 3), 3.814)
isequal(round(pvalue(lt), 4), 0.1485)

### X^2_S2
lt <- logrank_test(Surv(time, event) ~ group, data = leukemia,
                   type = "Tarone-Ware")
stopifnot(all(c(round(-statistic(lt, "linear")[1:2], 2),
                round(-statistic(lt, "linear")[3], 3)) == c(42.78, -26.42, -16.361)))
isequal(round(statistic(lt), 3), 4.104)
isequal(round(pvalue(lt), 4), 0.1285)
