
### compare output for examples from StatXact 6 manual with
### StatXact 6 output as given in the manual

library("coin")
set.seed(290875)
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)

### marginal homogenity

### StatXact 6 manual, page 315
presidents <- as.table(matrix(c(28, 7, 13, 27), nrow = 2,
    dimnames = list(BeforeTVDebate = c("Carter", "Reagan"),
                    AfterTVDebate = c("Carter", "Reagan"))))

bta <- mh_test(presidents)

# test statistic, page 306
stopifnot(isequal(round(sqrt(statistic(bta)), 3), 1.342))

# two-sided asymptotic p-value, page 306
stopifnot(isequal(round(pvalue(bta), 4), 0.1797))

# exact p-value, page 306
btMC <- mh_test(presidents, distribution = approximate(B = 10000))
pci <- attr(pvalue(btMC), "conf.int")
stopifnot(pci[1] < 0.2632 & pci[2] > 0.2632)


### StatXact 6 manual, page 308
endometrial_cancer <- as.table(matrix(c(6, 9, 9, 12,
                                        2, 4, 2,  1,
                                        3, 2, 3,  2,
                                        1, 1, 1, 1), nrow = 4,
                          dimnames = list(Cases = c(0, 0.2, 0.5125, 0.7),
                                          Controls = c(0, 0.2, 0.5125, 0.7))))

bta <- mh_test(endometrial_cancer,
                   scores = list(response = c(0, 0.2, 0.5125, 0.7)))

# test statistic, page 311
stopifnot(isequal(round(statistic(bta), 3), 3.735))

# two-sided asymptotic p-value, page 311
stopifnot(isequal(round(pvalue(bta), 4), 0.0002))

### StatXact 6 manual, page 312
pathologists <- as.table(matrix(c(22,  5,  0,  0, 0,
                                   2,  7,  2,  1, 0,
                                   2, 14, 36, 14, 3,
                                   0,  0,  0,  7, 0,
                                   0,  0,  0,  0, 3), nrow = 5,
    dimnames = list(Pathologist1 = paste("Level", 1:5, sep = "-"),
                    Pathologist2 = paste("Level", 1:5, sep = "-"))))

bta <- mh_test(pathologists,
                   scores = list(response = 1:5))

# test statistic, page 313
stopifnot(isequal(round(statistic(bta), 3), 1.152))

# two-sided asymptotic p-value, page 313
stopifnot(isequal(round(pvalue(bta), 4), 0.2492))

# exact p-value, page 313
btMC <- mh_test(pathologists, scores = list(response = 1:5),
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(btMC), "conf.int")
stopifnot(pci[1] < 0.3073 & pci[2] > 0.3073)


### Two independent samples

### StatXact 6 manual, page 340
load("bloodp.rda")

wta <- wilcox_test(bp ~ group, data = bloodp)

# test statistic, page 342
stopifnot(isequal(round(statistic(wta), 3), 1.720))

# two-sided asymptotic p-value, page 342
stopifnot(isequal(round(pvalue(wta), 4), 0.0853))

wte <- wilcox_test(bp ~ group, data = bloodp, distribution = "exact")

# two-sided exact p-value, page 342
stopifnot(isequal(round(pvalue(wte), 4), 0.0989))

wtel <- wilcox_test(bp ~ group, data = bloodp, distribution = "exact",
                     alternative = "greater")

# one-sided exact p-value, page 342
stopifnot(isequal(round(pvalue(wtel), 4), 0.0542))

# two-sided approximated p-value
wtMC <- wilcox_test(bp ~ group, data = bloodp,
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(wtMC), "conf.int")

stopifnot(pci[1] < pvalue(wte) & pci[2] > pvalue(wte))

# one-sided approximated p-value
wtMC <- wilcox_test(bp ~ group, data = bloodp, alternative = "greater",
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(wtMC), "conf.int")

stopifnot(pci[1] < pvalue(wtel) & pci[2] > pvalue(wtel))

# deal with ties as described in StatXact 6 manual, page 329ff
nta <- normal_test(bp ~ group, data = bloodp, ties.method = "average")

# test statistic, page 354
stopifnot(isequal(round(statistic(nta), 3), 1.789))

# two-sided asymptotic p-value, page 354
stopifnot(isequal(round(pvalue(nta), 4), 0.0737))

nte <- normal_test(bp ~ group, data = bloodp, distribution = "exact",
    ties.method = "average")

# two-sided exact p-value, page 354
stopifnot(isequal(round(pvalue(nte), 4), 0.0799))

ntel <- normal_test(bp ~ group, data = bloodp, distribution = "exact",
    ties.method = "average", alternative = "greater")

# one-sided exact p-value, page 354
stopifnot(isequal(round(pvalue(ntel), 4), 0.0462))

# two-sided approximated p-value
ntMC <- normal_test(bp ~ group, data = bloodp, ties.method = "average",
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(ntMC), "conf.int")

stopifnot(pci[1] < pvalue(nte) & pci[2] > pvalue(nte))

# one-sided approximated p-value
ntMC <- normal_test(bp ~ group, data = bloodp, alternative = "greater",
                    distribution = approximate(B = 10000),
                    ties.method = "average")
pci <- attr(pvalue(ntMC), "conf.int")

stopifnot(pci[1] < pvalue(ntel) & pci[2] > pvalue(ntel))

pta <- oneway_test(bp ~ group, data = bloodp)

# test statistic, page 394
stopifnot(isequal(round(statistic(pta), 3), 1.612))

# two-sided asymptotic p-value, page 394
stopifnot(isequal(round(pvalue(pta), 4), 0.1070))

pte <- oneway_test(bp ~ group, data = bloodp, distribution = "exact")

# two-sided exact p-value, page 394
stopifnot(isequal(round(pvalue(pte), 4), 0.1040))

ptel <- oneway_test(bp ~ group, data = bloodp, distribution = "exact",
                  alternative = "greater")

# one-sided exact p-value, page 394
stopifnot(isequal(round(pvalue(ptel), 4), 0.0564))

# two-sided approximated p-value
ptMC <- oneway_test(bp ~ group, data = bloodp,
                  distribution = approximate(B = 10000))
pci <- attr(pvalue(ptMC), "conf.int")

stopifnot(pci[1] < pvalue(pte) & pci[2] > pvalue(pte))

# one-sided approximated p-value
ptMC <- oneway_test(bp ~ group, data = bloodp, alternative = "greater",
                  distribution = approximate(B = 10000))
pci <- attr(pvalue(ptMC), "conf.int")

stopifnot(pci[1] < pvalue(ptel) & pci[2] > pvalue(ptel))


### StatXact 6 manual, page 345
load("employment.rda")

wta <- wilcox_test(Salary ~ Gender | Year, data = employment)

# test statistic, page 346
stopifnot(isequal(round(statistic(wta), 3), -1.897))

# two-sided asymptotic p-value, page 346
stopifnot(isequal(round(pvalue(wta), 4), 0.0578))

wte <- wilcox_test(Salary ~ Gender | Year, data = employment,
                   distribution = "exact")

# two-sided exact p-value, page 346
stopifnot(isequal(round(pvalue(wte), 4), 0.04))

wtel <- wilcox_test(Salary ~ Gender | Year, data = employment,
                    distribution = "exact", alternative = "less")

# one-sided exact p-value, page 346
stopifnot(isequal(round(pvalue(wtel), 4), 0.04))

# two-sided approximated p-value
wtMC <- wilcox_test(Salary ~ Gender | Year, data = employment,
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(wtMC), "conf.int")

stopifnot(pci[1] < pvalue(wte) & pci[2] > pvalue(wte))

# one-sided approximated p-value
wtMC <- wilcox_test(Salary ~ Gender | Year, data = employment,
                    alternative = "less",
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(wtMC), "conf.int")

stopifnot(pci[1] < pvalue(wtel) & pci[2] > pvalue(wtel))

nta <- normal_test(Salary ~ Gender | Year, data = employment,
    ties.method = "average")

# test statistic, page 358
stopifnot(isequal(round(statistic(nta), 3), -1.802))

# two-sided asymptotic p-value, page 358
stopifnot(isequal(round(pvalue(nta), 4), 0.0716))

# blocks not yet implemented
# nte <- normal_test(Salary ~ Gender | Year, data = employment,
#                   distribution = "exact")
#
# two-sided exact p-value, page 358
# stopifnot(isequal(round(pvalue(nte), 4), 0.04))
#
# ntel <- normal_test(Salary ~ Gender | Year, data = employment,
#                    distribution = "exact", alternative = "less")
#
# one-sided exact p-value, page 358
# stopifnot(isequal(round(pvalue(ntel), 4), 0.04))

# two-sided approximated p-value
ntMC <- normal_test(Salary ~ Gender | Year, data = employment,
    ties.method = "average", distribution = approximate(B = 10000))
pci <- attr(pvalue(ntMC), "conf.int")

stopifnot(pci[1] < 0.04 & pci[2] > 0.04)

# one-sided approximated p-value
ntMC <- normal_test(Salary ~ Gender | Year, data = employment,
                    alternative = "less", distribution = approximate(B = 10000),
                    ties.method = "average")
pci <- attr(pvalue(ntMC), "conf.int")

stopifnot(pci[1] < 0.04 & pci[2] > 0.04)

pta <- oneway_test(Salary ~ Gender | Year, data = employment)

# test statistic, page 396
stopifnot(isequal(round(statistic(pta), 3), -1.873))

# two-sided asymptotic p-value, page 396
stopifnot(isequal(round(pvalue(pta), 4), 0.0610))

pte <- oneway_test(Salary ~ Gender | Year, data = employment,
                 distribution = "exact")

# two-sided exact p-value, page 396
stopifnot(isequal(round(pvalue(pte), 4), 0.04))

ptel <- oneway_test(Salary ~ Gender | Year, data = employment,
                  distribution = "exact", alternative = "less")

# one-sided exact p-value, page 396
stopifnot(isequal(round(pvalue(ptel), 4), 0.04))

# two-sided approximated p-value
ptMC <- oneway_test(Salary ~ Gender | Year, data = employment,
                  distribution = approximate(B = 10000))
pci <- attr(pvalue(ptMC), "conf.int")

stopifnot(pci[1] < pvalue(pte) & pci[2] > pvalue(pte))

# one-sided approximated p-value
ptMC <- oneway_test(Salary ~ Gender | Year, data = employment,
                    alternative = "less",
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(ptMC), "conf.int")

stopifnot(pci[1] < pvalue(ptel) & pci[2] > pvalue(ptel))


### StatXact 9 manual, page 182
machines <- data.frame(cereal = c(10.8, 11.1, 10.4, 10.1, 11.3,
                                  10.8, 10.5, 11.0, 10.9, 10.8,
                                  10.7, 10.8),
                       machine = factor(rep(1:2, c(5, 7)),
                                        labels = c("Present", "New")))

ata <- ansari_test(cereal ~ machine, data = machines,
    ties.method = "average")

# test statistic, page 182
stopifnot(isequal(round(statistic(ata), 4), -1.9983))

# two-sided asymptotic p-value, page 182
stopifnot(isequal(round(pvalue(ata), 4), 0.0457))

ate <- ansari_test(cereal ~ machine, data = machines,
    ties.method = "average", distribution = "exact")

# two-sided exact p-value, page 182
stopifnot(isequal(round(pvalue(ate), 4), 0.0581))

atMC <- ansari_test(cereal ~ machine, data = machines,
    ties.method = "average", distribution = approximate(B = 10000))
pci <- attr(pvalue(atMC), "conf.int")

# two-sided approximated p-value, page 182
stopifnot(pci[1] < pvalue(ate) & pci[2] > pvalue(ate))

# Note: StatXact has '.LE.' here since *small* test statistics furnish evidence
#       for the alternative that sample 1 is *more* variable than sample 2, but
#       'coin' relieves the user from thinking about that detail

atag <- ansari_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater")

# one-sided asymptotic p-value, page 182
stopifnot(isequal(round(pvalue(atag), 4), 0.0228))

ateg <- ansari_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater",
    distribution = "exact")

# one-sided exact p-value, page 182
stopifnot(isequal(round(pvalue(ateg), 4), 0.0253))

atMC <- ansari_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater",
    distribution = approximate(B = 10000))
pci <- attr(pvalue(atMC), "conf.int")

# one-sided approximated p-value, page 182
stopifnot(pci[1] < pvalue(ateg) & pci[2] > pvalue(ateg))


### StatXact 9 manual, page 186

kta <- klotz_test(cereal ~ machine, data = machines,
    ties.method = "average")

# test statistic, page 186
stopifnot(isequal(round(statistic(kta), 4), 2.3082))

# two-sided asymptotic p-value, page 186
stopifnot(isequal(round(pvalue(kta), 4), 0.0210))

kte <- klotz_test(cereal ~ machine, data = machines,
    ties.method = "average", distribution = "exact")

# two-sided exact p-value, page 186
stopifnot(isequal(round(pvalue(kte), 4), 0.0101))

ktMC <- klotz_test(cereal ~ machine, data = machines,
    ties.method = "average", distribution = approximate(B = 10000))
pci <- attr(pvalue(ktMC), "conf.int")

# two-sided approximated p-value, page 186
stopifnot(pci[1] < pvalue(kte) & pci[2] > pvalue(kte))

ktag <- klotz_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater")

# one-sided asymptotic p-value, page 186
stopifnot(isequal(round(pvalue(ktag), 4), 0.0105))

kteg <- klotz_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater",
    distribution = "exact")

# one-sided exact p-value, page 186
stopifnot(isequal(round(pvalue(kteg), 4), 0.0101))

ktMC <- klotz_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater",
    distribution = approximate(B = 10000))
pci <- attr(pvalue(ktMC), "conf.int")

# one-sided approximated p-value, page 186
stopifnot(pci[1] < pvalue(kteg) & pci[2] > pvalue(kteg))


### StatXact 9 manual, page 190

mta <- mood_test(cereal ~ machine, data = machines,
    ties.method = "average")

# test statistic, page 190
stopifnot(isequal(round(statistic(mta), 4), 2.2715))

# two-sided asymptotic p-value, page 190
stopifnot(isequal(round(pvalue(mta), 4), 0.0231))

mte <- mood_test(cereal ~ machine, data = machines,
    ties.method = "average", distribution = "exact")

# two-sided exact p-value, page 190
stopifnot(isequal(round(pvalue(mte), 4), 0.0202))

mtMC <- mood_test(cereal ~ machine, data = machines,
    ties.method = "average", distribution = approximate(B = 10000))
pci <- attr(pvalue(mtMC), "conf.int")

# two-sided approximated p-value, page 190
stopifnot(pci[1] < pvalue(mte) & pci[2] > pvalue(mte))

mtag <- mood_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater")

# one-sided asymptotic p-value, page 190
stopifnot(isequal(round(pvalue(mtag), 4), 0.0116))

mteg <- mood_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater",
    distribution = "exact")

# one-sided exact p-value, page 190
stopifnot(isequal(round(pvalue(mteg), 4), 0.0126))

mtMC <- mood_test(cereal ~ machine, data = machines,
    ties.method = "average", alternative = "greater",
    distribution = approximate(B = 10000))
pci <- attr(pvalue(mtMC), "conf.int")

# one-sided approximated p-value, page 190
stopifnot(pci[1] < pvalue(mteg) & pci[2] > pvalue(mteg))


### StatXact 9 manual, page 195
load("FAILURE.rda")

cta <- conover_test(ftime ~ tyretype, data = failure)

# test statistic, page 195
stopifnot(isequal(round(statistic(cta), 4), -1.5274))

# two-sided asymptotic p-value, page 195
stopifnot(isequal(round(pvalue(cta), 4), 0.1267))

cte <- conover_test(ftime ~ tyretype, data = failure,
    distribution = "exact")

# two-sided exact p-value, page 195
stopifnot(isequal(round(pvalue(cte), 4), 0.1300))

ctMC <- conover_test(ftime ~ tyretype, data = failure,
    distribution = approximate(B = 10000))
pci <- attr(pvalue(ctMC), "conf.int")

# two-sided approximated p-value, page 195
stopifnot(pci[1] < pvalue(cte) & pci[2] > pvalue(cte))

ctag <- conover_test(ftime ~ tyretype, data = failure,
    alternative = "less")

# one-sided asymptotic p-value, page 195
stopifnot(isequal(round(pvalue(ctag), 4), 0.0633))

cteg <- conover_test(ftime ~ tyretype, data = failure,
    alternative = "less", distribution = "exact")

# one-sided exact p-value, page 195
stopifnot(isequal(round(pvalue(cteg), 4), 0.0634))

ctMC <- conover_test(ftime ~ tyretype, data = failure,
    alternative = "less", distribution = approximate(B = 10000))
pci <- attr(pvalue(ctMC), "conf.int")

# one-sided approximated p-value, page 195
stopifnot(pci[1] < pvalue(cteg) & pci[2] > pvalue(cteg))


### StatXact 6 manual, 413
load("lungcancer.rda")

### NOTE: StatXact uses another tie handling method, see Callaert (2003)
### and the examples below (at the end of this file)

lta <- logrank_test(Surv(time, cens) ~ group, data = lungcancer, ties = "aver")

# test statistic, page 415
isequal(round(statistic(lta), 3), 2.949) # no "-", StatXact 9, p. 215

# two-sided asymptotic p-value, page 415
isequal(round(pvalue(lta), 4), 0.0032)

lte <- logrank_test(Surv(time, cens) ~ group, data = lungcancer,
                    distribution = "exact")

# two-sided exact p-value, page 415
stopifnot(isequal(round(pvalue(lte), 4), 0.0010))

ltel <- logrank_test(Surv(time, cens) ~ group, data = lungcancer,
                     distribution = "exact", alternative = "great")

# one-sided exact p-value, page 415
stopifnot(isequal(round(pvalue(ltel), 4), 0.0010))

# two-sided approximated p-value
ltMC <- logrank_test(Surv(time, cens) ~ group, data = lungcancer,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(ltMC), "conf.int")

stopifnot(pci[1] < pvalue(lte) & pci[2] > pvalue(lte))

# one-sided approximated p-value
ltMC <- logrank_test(Surv(time, cens) ~ group, data = lungcancer,
                     alternative = "greater",
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(ltMC), "conf.int")

stopifnot(pci[1] < pvalue(ltel) & pci[2] > pvalue(ltel))


### StatXact manual 6, page 418
srv <- data.frame(time = c(3, 5, 7, 8, 18, 12, 19, 20, 20, 33),
                  event = c(0, 0, 1, 1, 1, 1, 1, 1, 0, 0),
                  treatment = gl(2, 5),
                  gender = factor(c("M", "M", "F", "F", "M",
                                    "M", "M", "F", "F", "F")))

# page 419
logrank_test(Surv(time, event) ~ treatment | gender, data = srv,
             ties = "average")

pvalue(logrank_test(Surv(time, event) ~ treatment | gender, data = srv,
                    ties = "average",
                    distribution = approximate(B = 10000)))


### StatXact 6 manual, page 448
hypnosis <- data.frame(skin = c(23, 58, 11, 24, 34,
                               23, 53, 10, 20, 40,
                               23, 54, 22, 21, 22),
    subject = factor(rep(c(1, 2, 3), rep(5, 3))),
    treatment = factor(rep(1:5, 3),
        labels = c("Fear", "Happiness", "Depression",
                   "Calmness", "Agitation")))

fta <- friedman_test(skin ~  treatment | subject, data = hypnosis)

# test statistic, page 449
stopifnot(isequal(round(statistic(fta), 3), 9.153))

# two-sided asymptotical p-value, page 449
stopifnot(isequal(round(pvalue(fta), 4), 0.0574))

# two-sided approximated p-value

ftMC <- friedman_test(skin ~  treatment | subject, data = hypnosis,
                      distribution = approximate(B = 10000))
pci <- attr(pvalue(ftMC), "conf.int")

stopifnot(pci[1] < 0.0268 & pci[2] > 0.0268)


### StatXact 6 manual, page 458

analgesic_eff <- data.frame(response = factor(
    c(0, 1, 1,
      0, 1, 1,
      1, 0, 1,
      0, 0, 0,
      0, 0, 1,
      0, 1, 1,
      1, 0, 1,
      0, 0, 1,
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      0, 0, 1), labels = c("success", "failure")),
    treatment = factor(rep(c("Placebo", "Aspirin", "NewDrug"), 12)),
    subject = factor(rep(1:12, rep(3, 12))))

bta <- mh_test(response ~  treatment | subject, data = analgesic_eff)

# asymptotic p-value, page 459 (no frame, see text!)
stopifnot(isequal(round(pvalue(bta), 3), 0.02))

# approximative p-value
btMC <- mh_test(response ~  treatment | subject,
    data = analgesic_eff, distribution = approximate(B = 10000))

pci <- attr(pvalue(btMC), "conf.int")
stopifnot(pci[1] < 0.026 & pci[2] > 0.026)


### StatXact 6 manual, page 467, Page test
cotton <- data.frame(strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                                  7.68, 7.57, 7.73, 8.15, 8.00,
                                  7.21, 7.80, 7.74, 7.87, 7.93),
                     potash = ordered(rep(c(144, 108, 72, 54, 36), 3),
                                      levels = c(144, 108, 72, 54, 36)),
                     block = factor(rep(1:3, rep(5, 3))))

fta <- friedman_test(strength ~ potash | block, data = cotton)

# test statistic, page 467
stopifnot(isequal(round(statistic(fta), 3), 2.656))

# asymptotical p-value, page 467
stopifnot(isequal(round(pvalue(fta), 4), 0.0079))

# approximate p-value
ftMC <- friedman_test(strength ~ potash | block, data = cotton,
                      distribution = approximate(B = 10000))

pci <- attr(pvalue(ftMC), "conf.int")
stopifnot(pci[1] < 0.005 & pci[2] > 0.005)

class(cotton$potash) <- "factor"

# approximate p-value
ftMC <- friedman_test(strength ~ potash | block, data = cotton,
                      distribution = approximate(B = 10000))

pci <- attr(pvalue(ftMC), "conf.int")
stopifnot(pci[1] < 0.0376 & pci[2] > 0.0376)


### StatXact 6 manual, page 486
tox <- data.frame(response = c(0,  1,  8, 10,
                               0,  0,  3,  3,  8,
                               5,  6,  7, 14, 14,
                               1,  1,  6,  7,  7, 7, 8, 8, 10,
                               7, 10, 11, 12, 13),
                  drug = factor(paste("Drug", c(rep(1, 4),
                                                rep(2, 5),
                                                rep(3, 5),
                                                rep(4, 9),
                                                rep(5, 5)))))


mta <- independence_test(response ~ drug, data = tox,
                         ytrafo = function(data)
                             trafo(data, numeric_trafo = median_trafo),
                         teststat = "quad")

# StatXact reports results of scaled Pearson statistic!
a <- factor(tox$response <= 7)
mta <-  chisq_test(a ~ drug, data = tox)

# test statistic, page 487
stopifnot(isequal(round(statistic(mta), 3), 4.317))

# asymptotic p-value, page 487
stopifnot(isequal(round(pvalue(mta), 4), 0.3648))

mtMC <- independence_test(response ~ drug, data = tox,
                          ytrafo = function(data)
                              trafo(data, numeric_trafo = median_trafo),
                          teststat = "quad",
                          distribution = approximate(B = 10000))

pci <- attr(pvalue(mtMC), "conf.int")
pvalue(mtMC)

stopifnot(pci[1] < 0.4289 & pci[2] > 0.4289)


kta <- kruskal_test(response ~ drug, data = tox)

# test statistic, page 493
stopifnot(isequal(round(statistic(kta), 3), 9.415))

# asymptotic p-value, page 493
stopifnot(isequal(round(pvalue(kta), 4), 0.0515))

pta <- oneway_test(response ~ drug, data = tox, teststat = "quad")

# test statistic, page 508
stopifnot(isequal(round(statistic(pta), 2), 10.98))

# asymptotic p-value, page 508
stopifnot(isequal(round(pvalue(pta), 4), 0.0268))


### StatXact 6 manual, page 509
brain <- data.frame(time = c(4,  5,  9, 12, 20, 25, 30,
                             1,  4,  9, 12, 15, 23, 30,
                             3,  7, 14, 20, 27, 30, 32, 50,
                             5, 15, 20, 31, 39, 47, 55, 67),
                    event = c(1, 1, 1, 1, 0, 1, 0,
                              1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 0, 0,
                              1, 1, 1, 1, 1, 1, 0,0),
                    treatment = factor(paste("trt",
                                             c(rep(1, 7), rep(2, 7),
                                               rep(3, 8), rep(3, 8)))))

lta <- logrank_test(Surv(time, event) ~ treatment, data = brain,
                    ties.method = "average-scores")

# test statistic, page 516
isequal(round(statistic(lta), 3), 5.012)

# asymptotic p-value, page 516
isequal(round(pvalue(lta), 4), 0.1709)


### StatXact 6 manual, page 519
oring <- data.frame(temp = c(66, 67, 67, 67, 68, 68, 70, 70, 72,
                             73, 75, 76, 76, 78, 79, 80, 81,
                             57, 58, 63, 70, 70,
                             75,
                             53),
                    incidents = ordered(c(rep("None", 17), rep("One", 5),
                                          rep("Two", 1), rep("Three", 1)),
                                        levels = c("None", "One", "Two", "Three")))

pta <- oneway_test(temp ~ incidents, data = oring)

# test statistic, page 524
stopifnot(isequal(round(statistic(pta), 3), -2.698))

stopifnot(isequal(round(pvalue(pta), 4), 0.007))


brain$treatment <- ordered(brain$treatment)

lta <- logrank_test(Surv(time, event) ~ treatment, data = brain,
                    ties.method = "average-scores")

# test statistic, page 536
isequal(round(statistic(lta), 3), 1.773)

# asymptotic p-value, page 536
isequal(round(pvalue(lta), 4), 0.0762)


### StatXact 6 manual, ...
a <- as.table(matrix(c(5, 5, 9, 1), nrow = 2,
    dimnames = list(response = c("Response", "No Response"),
                    drug = c("A", "B"))))

# exact p-value, page 671
stopifnot(isequal(round(fisher.test(a)$p.value, 4), 0.1409))

cta <- chisq_test(a)

# test statistic, page 673
stopifnot(isequal(round(statistic(cta), 2), 3.81))

# asymptotic p-value, page 673
stopifnot(isequal(round(pvalue(cta), 3), 0.051))

ctMC <- chisq_test(a, distribution = approximate(B = 10000))
pci <- attr(pvalue(ctMC), "conf.int")

stopifnot(pci[1] < 0.1409 & pci[2] > 0.1409)

ctMC <- cmh_test(a, distribution = approximate(B = 10000))
pci <- attr(pvalue(ctMC), "conf.int")

stopifnot(pci[1] < 0.1409 & pci[2] > 0.1409)


### StatXact 6 manual, page 793
csom <- as.table(matrix(c(17066, 48, 14464, 38, 788, 5, 126, 1, 37, 1), nrow = 2,
    dimnames = list(Malformation = c("Absent", "Present"),
                    Alcohol = c("0", "<1", "1-2", "3-5", ">=6"))))

lta <- lbl_test(csom, scores = list(Alcohol = 0:4))

# test statistic, page 796
stopifnot(isequal(round(statistic(lta), 3), -1.352))

# asymptotic p-value, page 796
stopifnot(isequal(round(pvalue(lta), 4), 0.1764))
stopifnot(isequal(round(pvalue(lbl_test(csom)), 4),
                  round(prop.trend.test(csom[2,], colSums(csom))$p.value, 4)))

ltMC <- lbl_test(csom, distribution = approximate(B = 100))

pci <- attr(pvalue(ltMC), "conf.int")
stopifnot(pci[1] < 0.179 & pci[2] > 0.179)

lta <- lbl_test(csom, scores = list(Alcohol = c(0, 0.5, 1.5, 4, 7)))

# test statistic, page 807
stopifnot(isequal(round(statistic(lta), 3), -2.563))

# asymptotic p-value, page 807
stopifnot(isequal(round(pvalue(lta), 4), 0.0104))
stopifnot(isequal(round(pvalue(lbl_test(csom,
                               scores = list(Alcohol = c(0, 0.5, 1.5, 4, 7)))), 4),
                  round(prop.trend.test(csom[2,], colSums(csom),
                        score = c(0, 0.5, 1.5, 4, 7))$p.value, 4)))


### StatXact manual 6, page 797
load("LANCET.rda")

lta <- lbl_test(lancet)

# test statistic, page 799
stopifnot(isequal(round(statistic(lta), 3), 4.335))

# asymptotic p-value, page 799
stopifnot(isequal(round(pvalue(lta), 4), 0))

ltMC <- lbl_test(lancet, distribution = approximate(B = 10000))

pci <- attr(pvalue(ltMC), "conf.int")
stopifnot(pci[1] < 0.0000045 & pci[2] > 0.0000045)


### StatXact manual, page 808
tumor <- data.frame(number = c( 0,  0,  0,  1,
                                7, 10,  6,  8,
                                0,  1,  0,  1,
                               11,  9, 13, 14,
                                1,  1,  1,  2,
                               29, 26, 28, 20),
                    tumor = factor(rep(c(rep("Present", 4), rep("Absent", 4)), 3),
                                   levels = c("Present", "Absent")),
                    dose = ordered(rep(paste(c(0, 1, 5, 50), "units"), 6),
                                   levels = paste(c(0, 1, 5, 50), "units")),
                    stratum = factor(rep(paste("Stratum", 3:5), rep(8, 3))))

lta <- lbl_test(tumor ~ dose | stratum, data = tumor, weights = ~ number,
    scores = list(dose = c(0, 1, 5, 50)))

# test statistic, page 812
stopifnot(isequal(round(statistic(lta), 3), 1.739))

# asymptotic p-value, page 812
stopifnot(isequal(round(pvalue(lta), 3), 0.082))

lta <- lbl_test(xtabs(number ~ dose + tumor + stratum, data = tumor),
                scores = list(dose = c(0, 1, 5, 50)))

# test statistic, page 812
stopifnot(isequal(round(statistic(lta), 3), 1.739))

# asymptotic p-value, page 812
stopifnot(isequal(round(pvalue(lta), 3), 0.082))

ltMC <- lbl_test(xtabs(number ~ dose + tumor + stratum, data = tumor),
                 scores = list(dose = c(0, 1, 5, 50)),
                 distribution = approximate(B = 10000))

pci <- attr(pvalue(ltMC), "conf.int")
stopifnot(pci[1] < 0.0769 & pci[2] > 0.0769)


### StatXact 6 manual, page 832
endo <- as.table(matrix(c(6, 9, 9, 12, 2, 4, 2, 1, 3, 2, 3, 2, 1, 1, 1, 1),
    nrow = 4, dimnames =
        list(Cases = c("0", "0.1-0.299", "0.3-0.625", ">0.625"),
             Controls = c("0", "0.1-0.299", "0.3-0.625", ">0.625"))))

bta <- mh_test(endo, scores = list(response = c(0, 0.2, 0.5125, 0.7)))

# test statistic, page 837
stopifnot(isequal(round(statistic(bta), 3), 3.735))

# asymptotic p-value, page 837
stopifnot(isequal(round(pvalue(bta), 4), 0.0002))

btMC <- mh_test(endo, scores = list(response = c(0, 0.2, 0.5125, 0.7)),
                    distribution = approximate(B = 10000))

print(pvalue(btMC))
pci <- attr(pvalue(btMC), "conf.int")
(pci[1] < 0.0001 & pci[2] > 0.0001)


### StatXact 6 manual, page 944
oral_lesions <- as.table(matrix(c(0, 8, 0, 0, 0, 0, 0, 1, 1,
                                  1, 1, 1, 1, 1, 1, 1, 0, 0,
                                  0, 8, 0, 0, 0, 0, 0, 1, 1), ncol = 3,
    dimnames = list(site = paste("site", 1:9),
                    region = c("Kerala", "Gujarat", "Andhra"))))

cta <- chisq_test(oral_lesions)

# test statistic, page 947
stopifnot(isequal(round(statistic(cta), 1), 22.1))

# asymptotic p-value, page 947
stopifnot(isequal(round(pvalue(cta), 2), 0.14))

# approximate p-value, page 947
ctMC <- cmh_test(oral_lesions, distribution = approximate(B = 10000))

pci <- attr(pvalue(ctMC), "conf.int")
stopifnot(pci[1] < 0.0269 & pci[2] > 0.0269)

### StatXact 6 manual, 984
dr <- as.table(matrix(c(100, 18, 50, 50, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1), nrow = 4,
    dimnames = list(dose = paste((1:4)*100, "mg"),
                    tox = c("mild", "moderate", "severe", "death"))))

lta <- lbl_test(dr, scores = list(dose = (1:4)*100, tox = 1:4))

# teststatistic, page 993
stopifnot(isequal(round(statistic(lta), 3), 1.807))

# asymptotical p-value, page 993
stopifnot(isequal(round(pvalue(lta), 4), 0.0708))

ltMC <- lbl_test(dr, distribution = approximate(B = 10000))
pci <- attr(pvalue(ltMC), "conf.int")

stopifnot(pci[1] < 0.0792 & pci[2] > 0.0792)


lta <- lbl_test(dr, scores = list(tox = c(1, 3, 9, 27)))

# teststatistic, page 993
stopifnot(isequal(round(statistic(lta), 3), 1.734))

# asymptotical p-value, page 993
stopifnot(isequal(round(pvalue(lta), 4), 0.0828))

ltMC <- lbl_test(dr, scores = list(tox = c(1, 3, 9, 27)),
                 distribution = approximate(B = 10000))

pci <- attr(pvalue(ltMC), "conf.int")
stopifnot(pci[1] < 0.078 & pci[2] > 0.078)


### StatXact 6 manual, page 1012
load("army.rda")

cta <- cmh_test(army)

# teststatistic, page 1014
stopifnot(isequal(round(statistic(cta), 2), 25.18))

# asymptotical p-value, page 1014
stopifnot(isequal(round(pvalue(cta), 6), 0.002774))


###
data(jobsatisfaction)

cta <- cmh_test(jobsatisfaction)

# teststatistic, page 1017
stopifnot(isequal(round(statistic(cta), 1), 10.2))

# asymptotical p-value, page 1017
stopifnot(isequal(round(pvalue(cta), 4), 0.3345))

# only Job.Satisfaction ordered
cta <- cmh_test(jobsatisfaction,
    scores = list(Job.Satisfaction = c(1, 3, 4, 5)))

# teststatistic, page 1018, is _wrong_ (but StatXact itself
# computes the correct result)
# (isequal(round(statistic(cta), 3), 9.226))
# Agresti, 2002, Table 7.12, page 297
stopifnot(isequal(round(statistic(cta), 4), 9.0342))

# asymptotical p-value, page 1018, is _wrong_ (but StatXact itself
# computes the correct result)
# (isequal(round(pvalue(cta), 4), 0.02643))
# Agresti, 2002, Table 7.12, page 297
stopifnot(isequal(round(pvalue(cta), 4), 0.0288))


lta <- lbl_test(jobsatisfaction,
    scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                  Income = c(3, 10, 20, 35)))

# teststatistic, page 1020
stopifnot(isequal(round(statistic(lta), 3), 2.481))

# asymptotical p-value, page 1020
stopifnot(isequal(round(pvalue(lta), 5), 0.01309))

lta <- cmh_test(jobsatisfaction,
    scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                  Income = c(3, 10, 20, 35)))

# teststatistic, page 1020
stopifnot(isequal(round(statistic(lta), 3), 2.481))

# asymptotical p-value, page 1020
stopifnot(isequal(round(pvalue(lta), 5), 0.01309))


### --------------------------------------------------------- ###

# some additional checks (always add new tests at the end because of the RNG's)

lta <- logrank_test(Surv(time, event) ~ treatment | gender, data = srv)
stopifnot(isequal(round(pvalue(lta), 4), 0.1224))

### example from Callaert (2003), AmStat 57, 214-217
exdata <- data.frame(time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
                     event = rep(TRUE, 15),
                     group = factor(c(rep(0, 7), rep(1, 8))))
p <- pvalue(logrank_test(Surv(time, event) ~ group, data = exdata,
                         distribution = exact()))
stopifnot(isequal(round(p, 4), 0.0505))
p <- pvalue(logrank_test(Surv(time, event) ~ group, data = exdata,
                         distribution = exact(), ties = "average"))
stopifnot(isequal(round(p, 4), 0.0468))


### symmetry problem, page 288
load("AIDS.rda")

tmp <- data.frame(y = c(AIDS$post, AIDS$pre),
                  x = gl(2, nrow(AIDS)),
                  block = factor(rep(1:nrow(AIDS), 2)))

wsa <- wilcoxsign_test(post ~ pre, data = AIDS, distribution = "asymptotic",
                       zero.method = "Wilcoxon")
wsa2 <- wilcoxsign_test(y ~ x | block, data = tmp, distribution = "asymptotic",
                        zero.method = "Wilcoxon")

stopifnot(all.equal(statistic(wsa), statistic(wsa2)))

# asymptotic p-value, page 290
stopifnot(isequal(round(statistic(wsa), 3), -2.896))

# asymptotic p-value, page 290
stopifnot(isequal(round(pvalue(wsa), 4), 0.0038))

wsa <- wilcoxsign_test(post ~ pre, data = AIDS, distribution = "asymptotic", alternative = "less",
                       zero.method = "Wilcoxon")

# asymptotic p-value, page 290
stopifnot(isequal(round(statistic(wsa), 3), -2.896))

# asymptotic p-value, page 290
stopifnot(isequal(round(pvalue(wsa), 4), 0.0019))

wse <- wilcoxsign_test(post ~ pre, data = AIDS, distribution = "exact",
                       zero.method = "Wilcoxon")

# exact p-value, page 290
stopifnot(isequal(round(pvalue(wse), 4), 0.0021))

wsa <- wilcoxsign_test(post ~ pre, data = AIDS, distribution = "exact", alternative = "less",
                       zero.method = "Wilcoxon")

# exact one-sided p-value, page 290
stopifnot(isequal(round(pvalue(wsa), 4), 0.0011))

### permutation test
diff <- AIDS$pre - AIDS$post
diff <- diff[abs(diff) > 0]
y <- as.vector(t(cbind(abs(diff) * (diff > 0), abs(diff) * (diff <= 0))))
block <- gl(length(diff), 2)
x <- factor(rep(c("pos", "neg"), length(diff)))

sta <- symmetry_test(y ~ x | block)

# asymptotic p-value, page 300
stopifnot(isequal(round(statistic(sta), 3), -1.707))
stopifnot(isequal(round(pvalue(sta), 4), 0.0878))

# approximate p-value, page 300
stMC <- symmetry_test(y ~ x | block, distribution = approximate(B = 10000))
pci <- attr(pvalue(stMC), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021)


### sign test
sa <- sign_test(post ~ pre, data = AIDS)
sa2 <- sign_test(y ~ x | block, data = tmp)

stopifnot(all.equal(statistic(sa), statistic(sa2)))

# asymptotic p-value, page 107 SX9 manual
stopifnot(isequal(round(statistic(sa), 3), -3))

# asymptotic p-value, page 107 SX9 manual
stopifnot(isequal(round(pvalue(sa), 4), 0.0027))

sa <- sign_test(post ~ pre, data = AIDS, alternative = "less")

# asymptotic p-value, page 107 SX9 manual
stopifnot(isequal(round(pvalue(sa), 4), 0.0013))

se <- sign_test(post ~ pre, data = AIDS, distribution = "exact")

# exact p-value, page 107 SX9 manual
stopifnot(isequal(round(pvalue(se), 4), 0.0042))

sa <- sign_test(post ~ pre, data = AIDS, distribution = "exact",
                alternative = "less")

# exact one-sided p-value, page 107 SX9 manual
stopifnot(isequal(round(pvalue(sa), 4), 0.0021))
