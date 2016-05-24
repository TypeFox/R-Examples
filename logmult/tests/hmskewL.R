# Artificial example to check the hmskewL model

library(logmult)
data(ocg1973)

tab <- array(ocg1973, dim=c(nrow(ocg1973), ncol(ocg1973), 2))

model <- hmskewL(tab[5:1, 5:1,], weighting="uniform", start=NA)
ass <- model$assoc

# First score for Farmers is slightly different from the original article
stopifnot(isTRUE(all.equal(round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,],
                                 matrix(c(-0.08, -0.2, -0.23, -0.11, 0.61,
                                           0.34,  0.3, -0.13, -0.51, 0), 5, 2),
                           check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(round(ass$row[,,1] * sqrt(ass$phi[2,1]), d=2)[5:1,],
                                 matrix(c(-0.08, -0.2, -0.23, -0.11, 0.61,
                                           0.34,  0.3, -0.13, -0.51, 0), 5, 2),
                           check.attributes=FALSE)))

model2 <- hmskewL(tab[5:1, 5:1,], weighting="uniform", layer.effect.skew="heterogeneous")
stopifnot(isTRUE(all.equal(model$assoc$phi, model2$assoc$phi)))
stopifnot(isTRUE(all.equal(model$assoc$row[,,1], model2$assoc$row[,,1])))


# EGP class of cohabiting spouses where one is 30-60 (last occupation for inactive persons)
# French Labour Force Surveys, 1969 and 2011
tab2 <- structure(c(261L, 43L, 21L, 5L, 7L, 26L, 5L, 16L, 17L, 7L, 1L,
                    483L, 394L, 215L, 53L, 58L, 117L, 37L, 185L, 232L, 104L, 11L,
                    565L, 457L, 528L, 139L, 116L, 201L, 19L, 469L, 788L, 368L, 14L,
                    148L, 195L, 399L, 306L, 96L, 213L, 50L, 321L, 1327L, 1344L, 123L,
                    17L, 9L, 11L, 3L, 33L, 37L, 4L, 12L, 13L, 11L, 1L, 165L, 70L,
                    105L, 77L, 573L, 878L, 55L, 78L, 188L, 181L, 11L, 9L, 3L, 34L,
                    7L, 23L, 36L, 1918L, 13L, 88L, 228L, 87L, 26L, 10L, 16L, 3L,
                    0L, 1L, 2L, 27L, 31L, 24L, 1L, 88L, 115L, 191L, 57L, 50L, 98L,
                    32L, 201L, 552L, 493L, 19L, 49L, 115L, 317L, 122L, 58L, 110L,
                    38L, 316L, 1301L, 1622L, 66L, 0L, 3L, 9L, 6L, 4L, 13L, 15L, 7L,
                    56L, 135L, 143L, 919L, 189L, 54L, 32L, 74L, 64L, 19L, 113L, 86L,
                    40L, 3L, 875L, 519L, 183L, 97L, 129L, 129L, 62L, 329L, 343L,
                    195L, 12L, 513L, 330L, 271L, 126L, 188L, 145L, 70L, 382L, 578L,
                    388L, 15L, 250L, 236L, 180L, 217L, 126L, 155L, 52L, 356L, 965L,
                    634L, 42L, 28L, 8L, 10L, 5L, 59L, 14L, 2L, 20L, 30L, 6L, 1L,
                    61L, 30L, 20L, 7L, 52L, 106L, 4L, 23L, 38L, 17L, 2L, 4L, 1L,
                    2L, 4L, 6L, 3L, 135L, 4L, 7L, 8L, 2L, 53L, 19L, 7L, 3L, 8L, 6L,
                    3L, 39L, 31L, 22L, 2L, 28L, 44L, 34L, 20L, 19L, 16L, 11L, 64L,
                    165L, 105L, 9L, 41L, 35L, 42L, 52L, 20L, 38L, 10L, 111L, 299L,
                    309L, 12L, 1L, 3L, 2L, 2L, 0L, 4L, 25L, 4L, 32L, 19L, 16L),
                  .Dim = c(11L, 11L, 2L), class="table",
                  .Dimnames = structure(list(H = c("I", "II", "IIIa", "IIIb", "IVa",
                                                   "IVb", "IVc", "V", "VI", "VIIa", "VIIb"),
                                             F = c("I", "II", "IIIa", "IIIb", "IVa", "IVb",
                                                   "IVc", "V", "VI", "VIIa", "VIIb"),
                                             T = c("1969", "2011")),
                  .Names = c("M", "W", "T")))

model2 <- hmskewL(tab2, start=NA)

stopifnot(isTRUE(all.equal(round(c(model2$assoc$phi), 2), c(0.18, 0.04, 0.18, 0.04))))
stopifnot(isTRUE(all.equal(round(c(model2$assoc$row), 2),
                           c(1.97, 1.38,  0.05, -0.24, -1.02, 0.57, -0.79, -0.16, -0.14, -1.32, -1.56,
                             0,   -0.03, -0.94,  0.98, -0.10, 0.71,  2.65, -0.86, -0.66, -0.77,  0.79))))

# Test anova
indep <- gnm(Freq ~ M*T + W*T, data=tab2, family=poisson)
anova(indep, model2, test="LR")
anova(indep, model2, test="Chisq")
