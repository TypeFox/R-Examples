### Regression tests for mid-p confidence intervals and mid-p-values

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)


###
### Mid-p confidence intervals
###


### Berry and Armitage (1995, p. 420)
mpci <- coin:::confint_midp(5, 20, 0.95)
stopifnot(isequal(round(mpci, 3), c(0.098, 0.470)))


### Agresti and Gottard (2001, p. 369)
mpci <- coin:::confint_midp(4, 10, 0.95)
stopifnot(isequal(round(mpci, 3), c(0.142, 0.709)))

mpci <- coin:::confint_midp(0, 10, 0.95)
stopifnot(isequal(round(mpci[1], 3), 0))

mpci <- coin:::confint_midp(10, 10, 0.95)
stopifnot(isequal(round(mpci[2], 3), 1))


### Newcombe (1998, p. 861, Tab. I)
mpci <- coin:::confint_midp(81, 263, 0.95)
stopifnot(isequal(round(mpci, 4), c(0.2544, 0.3658)))

mpci <- coin:::confint_midp(15, 148, 0.95)
stopifnot(isequal(round(mpci, 4), c(0.0601, 0.1581)))

mpci <- coin:::confint_midp(0, 20, 0.95)
stopifnot(isequal(round(mpci, 4), c(0.0000, 0.1391)))

mpci <- coin:::confint_midp(1, 29, 0.95)
stopifnot(isequal(round(mpci, 4), c(0.0017, 0.1585)))


###
### Mid-p-pvalues
###


### Data from Hwang and Yang (2001, p. 810)
tea <- matrix(c(3, 1,
                1, 3),
              nrow = 2, byrow = TRUE)

### Results from Hwang and Yang (2001, p. 810)
ct_e <- chisq_test(as.table(tea),
                   distribution = "exact")
stopifnot(isequal(round(pvalue(ct_e), 3), 0.486))
stopifnot(isequal(round(midpvalue(ct_e), 3), 0.257))

it_e_s <- independence_test(as.table(tea),
                            distribution = "exact",
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_e_s), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_s), midpvalue(ct_e)))

it_e_q <- independence_test(as.table(tea),
                            distribution = "exact",
                            teststat = "quad")
stopifnot(isequal(pvalue(it_e_q), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_q), midpvalue(ct_e)))

it_e_s_gr <- independence_test(as.table(tea),
                               distribution = "exact",
                               teststat = "scalar",
                               alternative = "greater")
stopifnot(isequal(round(pvalue(it_e_s_gr), 3), 0.243))
stopifnot(isequal(round(midpvalue(it_e_s_gr), 4), 0.1286))
                                 # p = 0.1285 according to Hwang and Yang (2001)

### Additional results: Monte Carlo
set.seed(290875)
ct_m <- chisq_test(as.table(tea),
                   distribution = "approximate")
(p <- pvalue(ct_m))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.486 & pci[2] > 0.486)
(mp <- midpvalue(ct_m))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.257 & mpci[2] > 0.257)

set.seed(290875)
it_m_s <- independence_test(as.table(tea),
                            distribution = approximate(B = 10000),
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_m_s), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_s), midpvalue(ct_m)))

set.seed(290875)
it_m_q <- independence_test(as.table(tea),
                            distribution = approximate(B = 10000),
                            teststat = "quad")
stopifnot(isequal(pvalue(it_m_q), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_q), midpvalue(ct_m)))

set.seed(290875)
it_m_s_gr <- independence_test(as.table(tea),
                               distribution = approximate(B = 10000),
                               teststat = "scalar",
                               alternative = "greater")
(p <- pvalue(it_m_s_gr))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.243 & pci[2] > 0.243)
(mp <- midpvalue(it_m_s_gr))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.1286 & mpci[2] > 0.1286)


### Data from Lydersen and Laake (2003, p. 3862)
davis <- matrix(c(3,  6,
                  2, 19),
                nrow = 2, byrow = TRUE)

### Results from Lydersen and Laake (2003, p. 3863, Tab. II)
ct_e <- chisq_test(as.table(davis),
                   distribution = "exact")
stopifnot(isequal(round(pvalue(ct_e), 4), 0.2860))
stopifnot(isequal(round(midpvalue(ct_e), 4), 0.1527))

it_e_s <- independence_test(as.table(davis),
                            distribution = "exact",
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_e_s), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_s), midpvalue(ct_e)))

it_e_q <- independence_test(as.table(davis),
                            distribution = "exact",
                            teststat = "quad")
stopifnot(isequal(pvalue(it_e_q), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_q), midpvalue(ct_e)))

### Additional results: Monte Carlo
set.seed(290875)
ct_m <- chisq_test(as.table(davis),
                   distribution = "approximate")
(p <- pvalue(ct_m))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.2860 & pci[2] > 0.2860)
(mp <- midpvalue(ct_m))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.1527 & mpci[2] > 0.1527)

set.seed(290875)
it_m_s <- independence_test(as.table(davis),
                            distribution = approximate(B = 10000),
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_m_s), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_s), midpvalue(ct_m)))

set.seed(290875)
it_m_q <- independence_test(as.table(davis),
                            distribution = approximate(B = 10000),
                            teststat = "quad")
stopifnot(isequal(pvalue(it_m_q), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_q), midpvalue(ct_m)))


### Data from Lydersen, Fagerland and Laake (2009, p. 1160, Tab. I)
cardiac <- matrix(c(1, 33,
                    7, 27),
                  nrow = 2, byrow = TRUE)

### Results from Lydersen, Fagerland and Laake (2009, pp. 1171--1172)
ct_e <- chisq_test(as.table(cardiac),
                   distribution = "exact")
stopifnot(isequal(round(pvalue(ct_e), 4), 0.0544))
stopifnot(isequal(round(midpvalue(ct_e), 4), 0.0297))

it_e_s <- independence_test(as.table(cardiac),
                            distribution = "exact",
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_e_s), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_s), midpvalue(ct_e)))

it_e_q <- independence_test(as.table(cardiac),
                            distribution = "exact",
                            teststat = "quad")
stopifnot(isequal(pvalue(it_e_q), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_q), midpvalue(ct_e)))

### Additional results: Monte Carlo
set.seed(290875)
ct_m <- chisq_test(as.table(cardiac),
                   distribution = approximate(B = 10000))
(p <- pvalue(ct_m))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.0544 & pci[2] > 0.0544)
(mp <- midpvalue(ct_m))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.0297 & mpci[2] > 0.0297)

set.seed(290875)
it_m_s <- independence_test(as.table(cardiac),
                            distribution = approximate(B = 10000),
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_m_s), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_s), midpvalue(ct_m)))

set.seed(290875)
it_m_q <- independence_test(as.table(cardiac),
                            distribution = approximate(B = 10000),
                            teststat = "quad")
stopifnot(isequal(pvalue(it_m_q), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_q), midpvalue(ct_m)))


### Data from Lydersen, Fagerland and Laake (2009, p. 1160, Tab. II)
exfoliative <- matrix(c( 0, 16,
                        15, 57),
                      nrow = 2, byrow = TRUE)

### Results from Lydersen, Fagerland and Laake (2009, p. 1173)
ct_e <- chisq_test(as.table(exfoliative),
                   distribution = "exact")
stopifnot(isequal(round(pvalue(ct_e), 4), 0.0629))
stopifnot(isequal(round(midpvalue(ct_e), 4), 0.0447))

it_e_s <- independence_test(as.table(exfoliative),
                            distribution = "exact",
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_e_s), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_s), midpvalue(ct_e)))

it_e_q <- independence_test(as.table(exfoliative),
                            distribution = "exact",
                            teststat = "quad")
stopifnot(isequal(pvalue(it_e_q), pvalue(ct_e)))
stopifnot(isequal(midpvalue(it_e_q), midpvalue(ct_e)))

### Additional results: Monte Carlo
set.seed(290875)
ct_m <- chisq_test(as.table(exfoliative),
                   distribution = approximate(B = 10000))
(p <- pvalue(ct_m))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.0629 & pci[2] > 0.0629)
(mp <- midpvalue(ct_m))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.0447 & mpci[2] > 0.0447)

set.seed(290875)
it_m_s <- independence_test(as.table(exfoliative),
                            distribution = approximate(B = 10000),
                            teststat = "scalar")
stopifnot(isequal(pvalue(it_m_s), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_s), midpvalue(ct_m)))

set.seed(290875)
it_m_q <- independence_test(as.table(exfoliative),
                            distribution = approximate(B = 10000),
                            teststat = "quad")
stopifnot(isequal(pvalue(it_m_q), pvalue(ct_m)))
stopifnot(isequal(midpvalue(it_m_q), midpvalue(ct_m)))


### Data from Fagerland, Lydersen and Laake (2013, p. 2, Tab. 1)
ahr <- matrix(c(1,  1,
                7, 12),
              nrow = 2, byrow = TRUE)

### Results from Fagerland, Lydersen and Laake (2013, p. 7, Tab. 6)
mt_e <- mh_test(as.table(ahr),
                distribution = "exact")
stopifnot(isequal(round(pvalue(mt_e), 4), 0.0703))
stopifnot(isequal(round(midpvalue(mt_e), 4), 0.0391))

st_e_s <- symmetry_test(as.table(ahr),
                        distribution = "exact",
                        teststat = "scalar")
stopifnot(isequal(pvalue(st_e_s), pvalue(mt_e)))
stopifnot(isequal(midpvalue(st_e_s), midpvalue(mt_e)))

st_e_q <- symmetry_test(as.table(ahr),
                        distribution = "exact",
                        teststat = "quad")
stopifnot(isequal(pvalue(st_e_q), pvalue(mt_e)))
stopifnot(isequal(midpvalue(st_e_q), midpvalue(mt_e)))

### Additional results: Monte Carlo
set.seed(290875)
mt_m <- mh_test(as.table(ahr),
                distribution = approximate(B = 10000))
(p <- pvalue(mt_m))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.0703 & pci[2] > 0.0703)
(mp <- midpvalue(mt_m))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.0391 & mpci[2] > 0.0391)

set.seed(290875)
st_m_s <- symmetry_test(as.table(ahr),
                        distribution = approximate(B = 10000),
                        teststat = "scalar")
stopifnot(isequal(pvalue(st_m_s), pvalue(mt_m)))
stopifnot(isequal(midpvalue(st_m_s), midpvalue(mt_m)))

set.seed(290875)
st_m_q <- symmetry_test(as.table(ahr),
                        distribution = approximate(B = 10000),
                        teststat = "quad")
stopifnot(isequal(pvalue(st_m_q), pvalue(mt_m)))
stopifnot(isequal(midpvalue(st_m_q), midpvalue(mt_m)))


### Data from Fagerland, Lydersen and Laake (2013, p. 2, Tab. 2)
therapy <- matrix(c(59,  6,
                    16, 80),
                  nrow = 2, byrow = TRUE)

### Results from Fagerland, Lydersen and Laake (2013, p. 7, Tab. 6)
mt_e <- mh_test(as.table(therapy),
                distribution = "exact")
stopifnot(isequal(round(pvalue(mt_e), 4), 0.0525))
stopifnot(isequal(round(midpvalue(mt_e), 4), 0.0347))

st_e_s <- symmetry_test(as.table(therapy),
                        distribution = "exact",
                        teststat = "scalar")
stopifnot(isequal(pvalue(st_e_s), pvalue(mt_e)))
stopifnot(isequal(midpvalue(st_e_s), midpvalue(mt_e)))

st_e_q <- symmetry_test(as.table(therapy),
                        distribution = "exact",
                        teststat = "quad")
stopifnot(isequal(pvalue(st_e_q), pvalue(mt_e)))
stopifnot(isequal(midpvalue(st_e_q), midpvalue(mt_e)))

### Additional results: Monte Carlo
set.seed(290875)
mt_m <- mh_test(as.table(therapy),
                distribution = approximate(B = 10000))
(p <- pvalue(mt_m))
pci <- attr(p, "conf.int")
stopifnot(pci[1] < 0.0525 & pci[2] > 0.0525)
(mp <- midpvalue(mt_m))
mpci <- attr(mp, "conf.int")
stopifnot(mpci[1] < 0.0347 & mpci[2] > 0.0347)

set.seed(290875)
st_m_s <- symmetry_test(as.table(therapy),
                        distribution = approximate(B = 10000),
                        teststat = "scalar")
stopifnot(isequal(pvalue(st_m_s), pvalue(mt_m)))
stopifnot(isequal(midpvalue(st_m_s), midpvalue(mt_m)))

set.seed(290875)
st_m_q <- symmetry_test(as.table(therapy),
                        distribution = approximate(B = 10000),
                        teststat = "quad")
stopifnot(isequal(pvalue(st_m_q), pvalue(mt_m)))
stopifnot(isequal(midpvalue(st_m_q), midpvalue(mt_m)))


### Data from Barnard (1989, p. 1470)
barnard50 <- matrix(c(5, 0,
                      0, 5),
                    nrow = 2, byrow = TRUE)

barnard51 <- matrix(c(5, 0,
                      1, 4),
                    nrow = 2, byrow = TRUE)

barnard52 <- matrix(c(5, 0,
                      2, 3),
                    nrow = 2, byrow = TRUE)

barnard53 <- matrix(c(5, 0,
                      3, 2),
                    nrow = 2, byrow = TRUE)

barnard54 <- matrix(c(5, 0,
                      4, 1),
                    nrow = 2, byrow = TRUE)

## barnard55 <- matrix(c(5, 0,
##                       5, 0),
##                     nrow = 2, byrow = TRUE)

barnard40 <- matrix(c(4, 1,
                      0, 5),
                    nrow = 2, byrow = TRUE)

barnard41 <- matrix(c(4, 1,
                      1, 4),
                    nrow = 2, byrow = TRUE)

barnard42 <- matrix(c(4, 1,
                      2, 3),
                    nrow = 2, byrow = TRUE)

barnard43 <- matrix(c(4, 1,
                      3, 2),
                    nrow = 2, byrow = TRUE)

barnard44 <- matrix(c(4, 1,
                      4, 1),
                    nrow = 2, byrow = TRUE)

barnard45 <- matrix(c(4, 1,
                      5, 0),
                    nrow = 2, byrow = TRUE)

barnard30 <- matrix(c(3, 2,
                      0, 5),
                    nrow = 2, byrow = TRUE)

barnard31 <- matrix(c(3, 2,
                      1, 4),
                    nrow = 2, byrow = TRUE)

barnard32 <- matrix(c(3, 2,
                      2, 3),
                    nrow = 2, byrow = TRUE)

barnard33 <- matrix(c(3, 2,
                      3, 2),
                    nrow = 2, byrow = TRUE)

barnard34 <- matrix(c(3, 2,
                      4, 1),
                    nrow = 2, byrow = TRUE)

barnard35 <- matrix(c(3, 2,
                      5, 0),
                    nrow = 2, byrow = TRUE)

barnard20 <- matrix(c(2, 3,
                      0, 5),
                    nrow = 2, byrow = TRUE)

barnard21 <- matrix(c(2, 3,
                      1, 4),
                    nrow = 2, byrow = TRUE)

barnard22 <- matrix(c(2, 3,
                      2, 3),
                    nrow = 2, byrow = TRUE)

barnard23 <- matrix(c(2, 3,
                      3, 2),
                    nrow = 2, byrow = TRUE)

barnard24 <- matrix(c(2, 3,
                      4, 1),
                    nrow = 2, byrow = TRUE)

barnard25 <- matrix(c(2, 3,
                      5, 0),
                    nrow = 2, byrow = TRUE)

barnard10 <- matrix(c(1, 4,
                      0, 5),
                    nrow = 2, byrow = TRUE)

barnard11 <- matrix(c(1, 4,
                      1, 4),
                    nrow = 2, byrow = TRUE)

barnard12 <- matrix(c(1, 4,
                      2, 3),
                    nrow = 2, byrow = TRUE)

barnard13 <- matrix(c(1, 4,
                      3, 2),
                    nrow = 2, byrow = TRUE)

barnard14 <- matrix(c(1, 4,
                      4, 1),
                    nrow = 2, byrow = TRUE)

barnard15 <- matrix(c(1, 4,
                      5, 0),
                    nrow = 2, byrow = TRUE)

## barnard00 <- matrix(c(0, 5,
##                       0, 5),
##                     nrow = 2, byrow = TRUE)

barnard01 <- matrix(c(0, 5,
                      1, 4),
                    nrow = 2, byrow = TRUE)

barnard02 <- matrix(c(0, 5,
                      2, 3),
                    nrow = 2, byrow = TRUE)

barnard03 <- matrix(c(0, 5,
                      3, 2),
                    nrow = 2, byrow = TRUE)

barnard04 <- matrix(c(0, 5,
                      4, 1),
                    nrow = 2, byrow = TRUE)

barnard05 <- matrix(c(0, 5,
                      5, 0),
                    nrow = 2, byrow = TRUE)

### Results from Barnard (1989, p. 1471, Tab. III; p. 1474, Tab. IV)
it50_e <- independence_test(as.table(barnard50),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it50_e), 4), 0.0040))
stopifnot(isequal(round(midpvalue(it50_e), 4), 0.0020))

it51_e <- independence_test(as.table(barnard51),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it51_e), 4), 0.0238))
stopifnot(isequal(round(midpvalue(it51_e), 4), 0.0119))

it52_e <- independence_test(as.table(barnard52),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it52_e), 4), 0.0833))
stopifnot(isequal(round(midpvalue(it52_e), 4), 0.0417))

it53_e <- independence_test(as.table(barnard53),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it53_e), 4), 0.2222))
stopifnot(isequal(round(midpvalue(it53_e), 4), 0.1111))

it54_e <- independence_test(as.table(barnard54),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it54_e), 4), 0.5000))
stopifnot(isequal(round(midpvalue(it54_e), 4), 0.2500))

## it55_e <- independence_test(as.table(barnard55),
##                             distribution = "exact",
##                             alternative = "greater")
## stopifnot(isequal(round(pvalue(it55_e), 4), 1.0000))
## stopifnot(isequal(round(midpvalue(it55_e), 4), 0.5000))

it40_e <- independence_test(as.table(barnard40),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it40_e), 4), 0.0238))
stopifnot(isequal(round(midpvalue(it40_e), 4), 0.0119))

it41_e <- independence_test(as.table(barnard41),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it41_e), 4), 0.1032))
stopifnot(isequal(round(midpvalue(it41_e), 4), 0.0536))

it42_e <- independence_test(as.table(barnard42),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it42_e), 4), 0.2619))
stopifnot(isequal(round(midpvalue(it42_e), 4), 0.1429))

it43_e <- independence_test(as.table(barnard43),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it43_e), 4), 0.5000))
stopifnot(isequal(round(midpvalue(it43_e), 4), 0.2917))

it44_e <- independence_test(as.table(barnard44),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it44_e), 4), 0.7778))
stopifnot(isequal(round(midpvalue(it44_e), 4), 0.5000))

it45_e <- independence_test(as.table(barnard45),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it45_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it45_e), 4), 0.7500))

it30_e <- independence_test(as.table(barnard30),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it30_e), 4), 0.0833))
stopifnot(isequal(round(midpvalue(it30_e), 4), 0.0417))

it31_e <- independence_test(as.table(barnard31),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it31_e), 4), 0.2619))
stopifnot(isequal(round(midpvalue(it31_e), 4), 0.1429))

it32_e <- independence_test(as.table(barnard32),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it32_e), 4), 0.5000))
stopifnot(isequal(round(midpvalue(it32_e), 4), 0.3016))

it33_e <- independence_test(as.table(barnard33),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it33_e), 4), 0.7381))
stopifnot(isequal(round(midpvalue(it33_e), 4), 0.5000))

it34_e <- independence_test(as.table(barnard34),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it34_e), 4), 0.9167))
stopifnot(isequal(round(midpvalue(it34_e), 4), 0.7083))

it35_e <- independence_test(as.table(barnard35),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it35_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it35_e), 4), 0.8889))

it20_e <- independence_test(as.table(barnard20),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it20_e), 4), 0.2222))
stopifnot(isequal(round(midpvalue(it20_e), 4), 0.1111))

it21_e <- independence_test(as.table(barnard21),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it21_e), 4), 0.5000))
stopifnot(isequal(round(midpvalue(it21_e), 4), 0.2917))

it22_e <- independence_test(as.table(barnard22),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it22_e), 4), 0.7381))
stopifnot(isequal(round(midpvalue(it22_e), 4), 0.5000))

it23_e <- independence_test(as.table(barnard23),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it23_e), 4), 0.8968))
stopifnot(isequal(round(midpvalue(it23_e), 4), 0.6984))

it24_e <- independence_test(as.table(barnard24),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it24_e), 4), 0.9762))
stopifnot(isequal(round(midpvalue(it24_e), 4), 0.8571))

it25_e <- independence_test(as.table(barnard25),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it25_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it25_e), 4), 0.9583))

it10_e <- independence_test(as.table(barnard10),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it10_e), 4), 0.5000))
stopifnot(isequal(round(midpvalue(it10_e), 4), 0.2500))

it11_e <- independence_test(as.table(barnard11),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it11_e), 4), 0.7778))
stopifnot(isequal(round(midpvalue(it11_e), 4), 0.5000))

it12_e <- independence_test(as.table(barnard12),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it12_e), 4), 0.9167))
stopifnot(isequal(round(midpvalue(it12_e), 4), 0.7083))

it13_e <- independence_test(as.table(barnard13),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it13_e), 4), 0.9762))
stopifnot(isequal(round(midpvalue(it13_e), 4), 0.8571))

it14_e <- independence_test(as.table(barnard14),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it14_e), 4), 0.9960))
stopifnot(isequal(round(midpvalue(it14_e), 4), 0.9464))

it15_e <- independence_test(as.table(barnard15),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it15_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it15_e), 4), 0.9881))

## it00_e <- independence_test(as.table(barnard00),
##                             distribution = "exact",
##                             alternative = "greater")
## stopifnot(isequal(round(pvalue(it00_e), 4), 1.0000))
## stopifnot(isequal(round(midpvalue(it00_e), 4), 0.5000))

it01_e <- independence_test(as.table(barnard01),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it01_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it01_e), 4), 0.7500))

it02_e <- independence_test(as.table(barnard02),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it02_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it02_e), 4), 0.8889))

it03_e <- independence_test(as.table(barnard03),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it03_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it03_e), 4), 0.9583))

it04_e <- independence_test(as.table(barnard04),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it04_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it04_e), 4), 0.9881))

it05_e <- independence_test(as.table(barnard05),
                            distribution = "exact",
                            alternative = "greater")
stopifnot(isequal(round(pvalue(it05_e), 4), 1.0000))
stopifnot(isequal(round(midpvalue(it05_e), 4), 0.9980))
