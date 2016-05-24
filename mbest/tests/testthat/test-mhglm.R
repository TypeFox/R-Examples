
context("mhglm")

library("lme4")

test_that("succeeds on sleepstudy", {
    model <- mhglm(Reaction ~ Days + (Days | Subject), data=sleepstudy)

    # fixef
    fixef0 <- c("(Intercept)" = 251.4, "Days" = 10.5)
    expect_that(round(fixef(model), 1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c(44.0, -1.4, -1.4, 2.3), 2, 2)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)", "Days")
    expect_that(round(vcov(model), 1), equals(vcov0))

    # VarCorr
    varcor0 <- matrix(c(565.5, 11.1, 11.1, 32.7), 2, 2)
    rownames(varcor0) <- colnames(varcor0) <- c("(Intercept)", "Days")
    varcor <- VarCorr(model)[["Subject"]]
    expect_that(attr(varcor, "stddev"), equals(sqrt(diag(varcor))))
    expect_that(as.vector(attr(varcor, "correlation")),
                equals(as.vector(t(varcor / attr(varcor, "stddev"))
                                 / attr(varcor, "stddev"))))
    attr(varcor, "stddev") <- NULL
    attr(varcor, "correlation") <- NULL
    expect_that(round(varcor, 1), equals(varcor0))

    # ranef
    ranef0 <- matrix(c( 2.8, -40.0, -38.4, 22.8, 21.5,  8.8,  16.4, -7.0,  -1.0,
                       34.7, -24.6, -12.3,  4.3, 20.6,  3.3, -24.7,  0.7,  12.1,
                        9.1,  -8.6,  -5.5, -4.7, -2.9, -0.2,  -0.2,  1.0, -10.6,
                        8.6,   1.1,   6.5, -3.0,  3.6,  0.9,   4.7, -1.0,   1.3),
                     18, 2)
    rownames(ranef0) <- as.character(c(308, 309, 310, 330, 331, 332, 333, 334,
                                       335, 337, 349, 350, 351, 352, 369, 370,
                                       371, 372))
    colnames(ranef0) <- c("(Intercept)", "Days")
    expect_that(round(as.matrix(ranef(model)[["Subject"]]), 1), equals(ranef0))
})

test_that("succeeds on sleepstudy in parallel", {
    
    model <- mhglm(Reaction ~ Days + (Days | Subject), data=sleepstudy,
                   control=list(parallel=TRUE))
    # fixef
    fixef0 <- c("(Intercept)" = 251.4, "Days" = 10.5)
    expect_that(round(fixef(model), 1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c(44.0, -1.4, -1.4, 2.3), 2, 2)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)", "Days")
    expect_that(round(vcov(model), 1), equals(vcov0))

    # VarCorr
    varcor0 <- matrix(c(565.5, 11.1, 11.1, 32.7), 2, 2)
    rownames(varcor0) <- colnames(varcor0) <- c("(Intercept)", "Days")
    varcor <- VarCorr(model)[["Subject"]]
    expect_that(attr(varcor, "stddev"), equals(sqrt(diag(varcor))))
    expect_that(as.vector(attr(varcor, "correlation")),
                equals(as.vector(t(varcor / attr(varcor, "stddev"))
                                 / attr(varcor, "stddev"))))
    attr(varcor, "stddev") <- NULL
    attr(varcor, "correlation") <- NULL
    expect_that(round(varcor, 1), equals(varcor0))

    # ranef
    ranef0 <- matrix(c( 2.8, -40.0, -38.4, 22.8, 21.5,  8.8,  16.4, -7.0,  -1.0,
                       34.7, -24.6, -12.3,  4.3, 20.6,  3.3, -24.7,  0.7,  12.1,
                        9.1,  -8.6,  -5.5, -4.7, -2.9, -0.2,  -0.2,  1.0, -10.6,
                        8.6,   1.1,   6.5, -3.0,  3.6,  0.9,   4.7, -1.0,   1.3),
                     18, 2)
    rownames(ranef0) <- as.character(c(308, 309, 310, 330, 331, 332, 333, 334,
                                       335, 337, 349, 350, 351, 352, 369, 370,
                                       371, 372))
    colnames(ranef0) <- c("(Intercept)", "Days")
    expect_that(round(as.matrix(ranef(model)[["Subject"]]), 1), equals(ranef0))
})


test_that("succeeds on cbpp", {
    suppressWarnings({
        model <- mhglm(cbind(incidence, size - incidence) ~ period + (period | herd),
                       data=cbpp, family=binomial)
    })

    # fixef
    fixef0 <- c("(Intercept)" = -1.2, "period2" = -0.8, "period3" = -0.9,
                "period4" = -1.1)
    expect_that(round(fixef(model), 1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c( 0.1, -0.1, -0.1, -0.1, -0.1, 0.2, 0.1, 0.0,
                      -0.1,  0.1,  0.2,  0.0, -0.1, 0.0, 0.0, 0.1), 4, 4)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)", "period2",
                                            "period3", "period4")
    expect_that(round(vcov(model), 1), equals(vcov0))

    # VarCorr
    varcor0 <- matrix(c( 0.7, -0.7, -0.4, -0.5, -0.7, 1.0, 1.0, 0.3,
                        -0.4,  1.0,  1.4,  0.0, -0.5, 0.3, 0.0, 0.3), 4, 4)
    rownames(varcor0) <- colnames(varcor0) <- c("(Intercept)", "period2",
                                                "period3", "period4")
    varcor <- VarCorr(model)[["herd"]]
    expect_that(attr(varcor, "stddev"), equals(sqrt(diag(varcor))))
    expect_that(as.vector(attr(varcor, "correlation")),
                equals(as.vector(t(varcor / attr(varcor, "stddev"))
                                 / attr(varcor, "stddev"))))
    attr(varcor, "stddev") <- NULL
    attr(varcor, "correlation") <- NULL
    expect_that(round(varcor, 1), equals(varcor0))

    # ranef
    ranef0 <- matrix(c(-0.1, -0.4,  0.5,  0.1,  0.3, -0.2,  0.9,  0.5,  0.0, -0.7,
                       -0.7,  0.0, -0.8,  1.2, -0.7,
                        0.8,  0.2, -0.6,  0.2, -0.8, -0.1, -0.6, -0.5, -0.2,  0.3,
                        1.0,  0.1,  0.7, -1.1,  0.7,
                        1.5, -0.1, -0.6,  0.5, -1.2, -0.5,  0.2, -0.3, -0.3, -0.2,
                        1.1,  0.1,  0.2, -0.7,  0.5,
                       -0.3,  0.3, -0.3, -0.2,  0.1,  0.3, -0.8, -0.3,  0.1,  0.6,
                        0.2,  0.0,  0.6, -0.7,  0.4),
                     15, 4)
    rownames(ranef0) <- as.character(1:15)
    colnames(ranef0) <- c("(Intercept)", "period2", "period3", "period4")
    expect_that(round(as.matrix(ranef(model)[["herd"]]), 1), equals(ranef0))
})

test_that("succeeds on cbpp in parallel", {
    suppressWarnings({
        model <- mhglm(cbind(incidence, size - incidence) ~ period + (period | herd),
                       data=cbpp, family=binomial,
                       control=list(parallel=TRUE))
    })

    # fixef
    fixef0 <- c("(Intercept)" = -1.2, "period2" = -0.8, "period3" = -0.9,
                "period4" = -1.1)
    expect_that(round(fixef(model), 1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c( 0.1, -0.1, -0.1, -0.1, -0.1, 0.2, 0.1, 0.0,
                      -0.1,  0.1,  0.2,  0.0, -0.1, 0.0, 0.0, 0.1), 4, 4)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)", "period2",
                                            "period3", "period4")
    expect_that(round(vcov(model), 1), equals(vcov0))

    # VarCorr
    varcor0 <- matrix(c( 0.7, -0.7, -0.4, -0.5, -0.7, 1.0, 1.0, 0.3,
                        -0.4,  1.0,  1.4,  0.0, -0.5, 0.3, 0.0, 0.3), 4, 4)
    rownames(varcor0) <- colnames(varcor0) <- c("(Intercept)", "period2",
                                                "period3", "period4")
    varcor <- VarCorr(model)[["herd"]]
    expect_that(attr(varcor, "stddev"), equals(sqrt(diag(varcor))))
    expect_that(as.vector(attr(varcor, "correlation")),
                equals(as.vector(t(varcor / attr(varcor, "stddev"))
                                 / attr(varcor, "stddev"))))
    attr(varcor, "stddev") <- NULL
    attr(varcor, "correlation") <- NULL
    expect_that(round(varcor, 1), equals(varcor0))

    # ranef
    ranef0 <- matrix(c(-0.1, -0.4,  0.5,  0.1,  0.3, -0.2,  0.9,  0.5,  0.0, -0.7,
                       -0.7,  0.0, -0.8,  1.2, -0.7,
                        0.8,  0.2, -0.6,  0.2, -0.8, -0.1, -0.6, -0.5, -0.2,  0.3,
                        1.0,  0.1,  0.7, -1.1,  0.7,
                        1.5, -0.1, -0.6,  0.5, -1.2, -0.5,  0.2, -0.3, -0.3, -0.2,
                        1.1,  0.1,  0.2, -0.7,  0.5,
                       -0.3,  0.3, -0.3, -0.2,  0.1,  0.3, -0.8, -0.3,  0.1,  0.6,
                        0.2,  0.0,  0.6, -0.7,  0.4),
                     15, 4)
    rownames(ranef0) <- as.character(1:15)
    colnames(ranef0) <- c("(Intercept)", "period2", "period3", "period4")
    expect_that(round(as.matrix(ranef(model)[["herd"]]), 1), equals(ranef0))
})

