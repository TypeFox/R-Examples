#############################################################################
# test-style_apa.R
# 
# Testing APA style functions
#############################################################################

# style.apa.t.test
#############################################################################

test_that("apa:t.test",
{
    t1 <- list(t.test(1:10))
    t2 <- list(t.test(1:10), 0.7819)
    t3 <- list(t.test(1:10, 2:11))

    expect_identical(style.apa.t.test(t1),
                     c("M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=5.50",
                       "t(9)=5.74",
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.t.test(t2),
                     c("M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=5.50",
                       "t(9)=5.74",
                       "p\\ifmmode<\\else\\textless\\fi.001",
                       "d=0.78"))
    expect_identical(style.apa.t.test(t3,
                                     print.estimate = TRUE,
                                     estimate.names = names(t3[[1]]$estimate)),
                     c("M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=5.50",
                       "M\\ifmmode_{y}\\else\\textsubscript{y}\\fi=6.50",
                       "t(18)=-0.74",
                       "p=.470"))
    expect_identical(style.apa.t.test(t3,
                                     print.estimate = TRUE,
                                     estimate.names = c("Mx", "My")),
                     c("Mx=5.50",
                       "My=6.50",
                       "t(18)=-0.74",
                       "p=.470"))
})

# style.apa.fisher
#############################################################################

test_that("apa:Fisher exact test",
{
    # example copied from help file
    Job <- matrix(c(1,2,1,0, 3,3,6,1,10,10,14,9, 6,7,12,11), 
                  4, 4, 
                  dimnames = list(income = c("< 15k", 
                                             "15-25k", 
                                             "25-40k", 
                                             "> 40k"), 
                                  satisfaction = c("VeryD", 
                                                   "LittleD",
                                                   "ModerateS", 
                                                   "VeryS")))
    Convictions <- matrix(c(2, 10, 15, 3), 
                          nrow = 2, 
                          dimnames = list(c("Dizygotic", 
                                            "Monozygotic"),
                                          c("Convicted", 
                                            "Not convicted")))

    t1 <- list(fisher.test(Job))
    t2 <- list(fisher.test(Convictions, alternative = "less"))

    expect_identical(style.apa.fisher(t1),
                     "p=.783")
    expect_identical(style.apa.fisher(t2),
                     c("odds ratio=0.05", 
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.fisher(t2,
                                      print.estimate = TRUE,
                                      estimate.names = "oddsratio"),
                     c("oddsratio=0.05", 
                       "p\\ifmmode<\\else\\textless\\fi.001"))
})

# style.apa.pearson
#############################################################################

test_that("apa: Pearson's product-moment correlation",
{
    # example copied from help file
    x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
    y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)

    t1 <- list(cor.test(x, y))

    expect_identical(style.apa.pearson(t1),
                     c("r=0.57", 
                       "t(7)=1.84", 
                       "p=.108"))
})

# style.apa.chisquared
#############################################################################

test_that("apa: chi-squared test",
{
    # example copied from help file
    M <- as.table(rbind(c(762, 327, 468), 
                        c(484, 239, 477)))
    dimnames(M) <- list(gender = c("M",
                                   "F"), 
                        party = c("Democrat",
                                  "Independent", 
                                  "Republican"))

    t1 <- list(chisq.test(M))

    expect_identical(style.apa.chisquared(t1),
                     c("\\ifmmode\\chi\\else\\(\\chi\\)\\fi\\ifmmode^{2}\\else\\textsuperscript{2}\\fi(2, N=2757)=30.07",
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.chisquared(t1,
                                          chi.n.name = 'n'),
                     c("\\ifmmode\\chi\\else\\(\\chi\\)\\fi\\ifmmode^{2}\\else\\textsuperscript{2}\\fi(2, n=2757)=30.07",
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.chisquared(t1,
                                          chi.n.name = NULL),
                     c("\\ifmmode\\chi\\else\\(\\chi\\)\\fi\\ifmmode^{2}\\else\\textsuperscript{2}\\fi(2, 2757)=30.07",
                       "p\\ifmmode<\\else\\textless\\fi.001"))
})

# style.apa.shapiro
#############################################################################

test_that("apa: shapiro test",
{
    set.seed(50)
    t1 <- list(shapiro.test(rnorm(100, mean = 5, sd = 3)))

    expect_identical(style.apa.shapiro(t1),
                     c("W=0.99",
                       "p=.932"))
})

# style.apa.bartlett
#############################################################################

test_that("apa: bartlett test",
{
    t1 <- list(bartlett.test(count ~ spray, data = InsectSprays))

    expect_identical(style.apa.bartlett(t1),
                     c("\\ifmmode\\chi\\else\\(\\chi\\)\\fi\\ifmmode^{2}\\else\\textsuperscript{2}\\fi(5)=25.96",
                       "p\\ifmmode<\\else\\textless\\fi.001"))
})

# style.apa.friedman
#############################################################################

test_that("apa: friedman test",
{
    wb <- aggregate(warpbreaks$breaks, 
                    by = list(w = warpbreaks$wool, 
                              t = warpbreaks$tension), 
                    FUN = mean)
    t1 <- list(friedman.test(x ~ w | t, data = wb))

    expect_identical(style.apa.bartlett(t1),
                     c("\\ifmmode\\chi\\else\\(\\chi\\)\\fi\\ifmmode^{2}\\else\\textsuperscript{2}\\fi(1)=0.33",
                       "p=.564"))
})

# style.apa.ks.test
#############################################################################

test_that("apa: Kolmogorov-Smirnow test",
{
    set.seed(50)
    t1 <- list(ks.test(rnorm(50), runif(30), alternative = "l"))

    expect_identical(style.apa.ks.test(t1),
                     c("D\\ifmmode^{-}\\else\\textsuperscript{-}\\fi=0.06",
                       "p=.874"))
})

# style.apa.summary.aovlist
#############################################################################

test_that("apa: one way anova",
{
    t1 <- list(summary(aov(count ~ spray, data = InsectSprays)))

    expect_identical(style.apa.summary.aovlist(t1),
                     c("F(5, 66)=34.70",
                       "p\\ifmmode<\\else\\textless\\fi.001"))
})

test_that("apa: two way anova with interaction",
{
    set.seed(50)
    myd <- data.frame(value = c(rnorm(20), rnorm(20) + 1),
                      effect1 = factor(rep(c(1, 2), times = 20)),
                      effect2 = factor(rep(c(1, 2), each = 20)))

    t1 <- list(summary(aov(value ~ effect1 * effect2, data = myd)))

    expect_identical(style.apa.summary.aovlist(t1, effect = 1),
                     c("F(1, 36)=1.00",
                       "p=.324"))
    expect_identical(style.apa.summary.aovlist(t1, effect = 2),
                     c("F(1, 36)=27.19", 
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.summary.aovlist(t1, effect = 3),
                     c("F(1, 36)=2.43", 
                       "p=.128"))
})

# style.apa.summary.lm
#############################################################################

test_that("apa: regression analysis -> summary.lm",
{
    set.seed(50)

    ## regression analysis
    response <- rnorm(100)
    predictor1 <- rnorm(100) * 3/4 + response * 1/4
    predictor2 <- rnorm(100) * 1/2 + response * 1/4 + predictor1 * 1/4
    
    myd <- data.frame(response,
                      predictor1,
                      predictor2)
    
    fit <- lm(response ~ predictor1 + predictor2, data = myd)
    fitsm <- list(summary(fit))
    
    expect_identical(style.apa.summary.lm(fitsm),
                     c("R\\ifmmode^{2}\\else\\textsuperscript{2}\\fi=0.20", 
                       "F(2, 97)=12.40", 
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.summary.lm(fitsm, result = "model", r.squared = "adjusted"),
                     c("R\\ifmmode^{2}\\else\\textsuperscript{2}\\fi=0.19", 
                       "F(2, 97)=12.40", 
                       "p\\ifmmode<\\else\\textless\\fi.001"))
    expect_identical(style.apa.summary.lm(fitsm, result = 1),
                     c("b=-0.11", 
                       "t(97)=-1.22", 
                       "p=.225"))
    expect_identical(style.apa.summary.lm(fitsm, result = 2),
                     c("b=0.19", 
                       "t(97)=1.53", 
                       "p=.129"))
    expect_identical(style.apa.summary.lm(fitsm, result = "equation"),
                     "\\ifmmode\\hat{y}\\else y\\textsuperscript{\\textasciicircum}\\fi = -0.11+0.19*predictor1+0.65*predictor2")
})
