library("papeR")
context("prettify functions")

set.seed(1234)

################################################################################
## Test computation of CIs when data is part of the call
## (i.e. not only a link is passed but really the data)
################################################################################

test_that("computation of CIs when data is part of the call", {
    model_fit <- function(my_data, model_class) {
        do.call(model_class, list(y ~ x, data = my_data))
    }

    for (model_class in c("lm", "glm")) {
        x <- rnorm(100)
        y <- rnorm(100, mean = 2 * x)
        data <- data.frame(y = y, x = x)

        ## fit model with data argument
        mod <- do.call(model_class, list(y ~ x, data = data))
        psm1 <- prettify(summary(mod))
        rm(data)
        psm1a <- prettify(summary(mod))

        ## fit model without data argument
        mod2 <- do.call(model_class, list(y ~ x))
        psm2 <- prettify(summary(mod2))

        ## fit model in different environment
        mod3 <- model_fit(data.frame(y = y, x = x), model_class)
        psm3 <- prettify(summary(mod3))

        ## change data and compute summary
        x <- rnorm(100)
        y <- rnorm(100, mean = 2 * x)
        data <- data.frame(y = y, x = x)

        psm4 <- prettify(summary(mod))

        expect_equal(psm1, psm1a)
        expect_equal(psm1, psm2)
        expect_equal(psm1, psm3)
        expect_equal(psm1, psm4)
    }
})

################################################################################
## Test computation of CIs when data is *not* part of the call
################################################################################

rm(list = ls())
model_fit2 <- function(my_data, model_class) {
    switch(model_class,
           lm = lm(y ~ x, data = my_data),
           glm = glm(y ~ x, data = my_data),
           coxph = coxph(Surv(y, event) ~ x, data = my_data),
           lme = lme(y ~ x, random = ~ 1 | group, data = my_data),
           lmer = lmer(y ~ x + (1 | group), data = my_data))
}

fit_model <- function(model_class =  c("lm", "glm", "coxph", "lme", "lmer")) {
    x <- rnorm(100)
    y <- rnorm(100, mean = 2 * x)
    data <- data.frame(y = y, x = x)

    if (model_class %in% c("lme", "lmer")) {
        group <- sample(1:2, 100, replace = TRUE)
        data$group <- group
    }

    if (model_class %in% "coxph") {
        event <- as.logical(sample(0:1, 100, replace = TRUE))
        data$event <- event
        data$y <- exp(y)
    }

    ## fit model with data argument
    mod <- switch(model_class,
                  lm = lm(y ~ x, data = data),
                  glm = glm(y ~ x, data = data),
                  coxph = coxph(Surv(y, event) ~ x, data = data),
                  lme = lme(y ~ x, random = ~ 1 | group, data = data),
                  lmer = lmer(y ~ x + (1 | group), data = data))
    psm1 <- prettify(summary(mod))
    data_tmp <- data
    rm(data)
    psm1a <- try(prettify(summary(mod)), silent = TRUE)

    ## fit model without data argument
    mod2 <- switch(model_class,
                   lm = lm(y ~ x),
                   glm = glm(y ~ x),
                   coxph = coxph(Surv(y, event) ~ x),
                   lme = lme(y ~ x, random = ~ 1 | group),
                   lmer = lmer(y ~ x + (1 | group)))
    psm2 <- prettify(summary(mod2))

    ## fit model in different environment
    mod3 <- model_fit2(data_tmp, model_class)
    psm3 <- try(prettify(summary(mod3)), silent = TRUE)

    ## change data and compute summary
    x <- rnorm(100)
    y <- rnorm(100, mean = 2 * x)
    data <- data.frame(y = y, x = x)

    if (model_class %in% c("lme", "lmer")) {
        group <- sample(1:2, 100, replace = TRUE)
        data$group <- group
    }

    if (model_class %in% "coxph") {
        event <- as.logical(sample(0:1, 100, replace = TRUE))
        data$event <- event
        data$y <- exp(y)
    }

    psm4 <- prettify(summary(mod))

    ret <- list(psm1, psm1a, psm2, psm3, psm4)
    names(ret) <- c("standard", "data removed from global environment",
                    "no data argument in call", "local environment",
                    "data in global environment changed")
    return(ret)
}


## check lm interface
test_that("lm interface works", {
    expect_warning(res <- fit_model("lm"))
    expect_equivalent(res[[1]], res[[3]])
    expect_equivalent(res[[1]], res[[4]])
    ## differences in CIs as different data is used
    expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[1], "lower")
    expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[2], "upper")
    ## CI dropped. Other values equal
    expect_equivalent(res[[1]][, -(3:4)], res[[2]])
})

## check glm interface
test_that("glm interface works", {
    expect_warning(res <- fit_model("glm"))
    expect_equivalent(res[[1]], res[[3]])
    expect_equivalent(res[[1]], res[[4]])
    ## differences in CIs as different data is used
    expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[1], "lower")
    expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[2], "upper")
    ## CI dropped. Other values equal
    expect_equivalent(res[[1]][, -(3:4)], res[[2]])
})

### check lme interface
if (require("nlme")) {
    test_that("nlme interface works", {
        expect_warning(res <- fit_model("lme"))
        expect_equivalent(res[[1]], res[[3]])
        expect_equivalent(res[[1]], res[[4]])
        ## differences in CIs as different data is used
        expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[1], "lower")
        expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[2], "upper")
        ## CI dropped. Other values equal
        stopifnot(all.equal(res[[1]][, -(3:4)], res[[2]], check.attributes = FALSE))
    })
}

## check coxph interfaces
if (require("survival")) {
    test_that("survival interface works", {
        expect_warning(res <- fit_model("coxph"))
        expect_equivalent(res[[1]], res[[3]])
        expect_is(res[[2]], "try-error")
        expect_is(res[[4]], "try-error")
        ## differences in CIs as different data is used
        expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[1], "lower")
        expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[2], "upper")
    })
}

## check lmer interface
if (require("lme4")) {
    test_that("lme4 interface works", {
        expect_warning(res <- fit_model("lmer"))
        if (packageDescription("lme4")$Version >= 1) {
            expect_equivalent(res[[1]], res[[3]])
            expect_is(res[[2]], "try-error")
            expect_is(res[[4]], "try-error")
            ## differences in CIs as different data is used
            expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[1], "lower")
            expect_match(all.equal(res[[1]], res[[5]], check.attributes = FALSE)[2], "upper")
        }
    })
}

################################################################################
## Test OR
################################################################################

test_that("OR are included", {
    x <- rnorm(100)
    y <- rbinom(100, 1, make.link("logit")$linkinv(x * 2))
    data <- data.frame(x, y)
    mod <- glm(y ~ x, data = data, family = binomial)
    ps <- prettify(summary(mod))
    expect_true("Odds Ratio" %in% names(ps))
})

################################################################################
## Test specification of CIs
################################################################################

fit_model <- function(model_class =  c("lm", "glm", "coxph", "lme", "lmer")) {
    x <- rnorm(100)
    y <- rnorm(100, mean = 2 * x)
    data <- data.frame(y = y, x = x)

    if (model_class %in% c("lme", "lmer")) {
        group <- sample(1:2, 100, replace = TRUE)
        data$group <- group
    }

    if (model_class %in% "coxph") {
        event <- as.logical(sample(0:1, 100, replace = TRUE))
        data$event <- event
        data$y <- exp(y)
    }

    ## fit model with data argument
    mod <- switch(model_class,
                  lm = lm(y ~ x, data = data),
                  glm = glm(y ~ x, data = data),
                  coxph = coxph(Surv(y, event) ~ x, data = data),
                  lme = lme(y ~ x, random = ~ 1 | group, data = data),
                  lmer = lmer(y ~ x + (1 | group), data = data))
    return(list(data = data, model = mod))
}

if (require("nlme") && require("lme4") && packageDescription("lme4")$Version >= 1) {
    test_that("CIs can be hand specified", {
        for (model_class in c("lm", "glm", "lme", "lmer")) {
            ## fit model
            RES <- fit_model(model_class)
            mod <- RES$mod
            data <- RES$data
            CI <- matrix(c(1, 2, 3, 4), ncol = 2)
            ps <- prettify(summary(mod), confint = CI)
            expect_equivalent(as.matrix(ps[,3:4]), CI, info = model_class)
        }
    })
}

################################################################################
## Test anova
################################################################################

if (require("nlme") && require("survival") && require("lme4") && packageDescription("lme4")$Version >= 1) {
    test_that("prettify.anova works", {
        for (model_class in c("lm", "glm", "lme", "lmer", "coxph")) {
            ## fit model
            RES <- fit_model(model_class)
            mod <- RES$mod
            data <- RES$data
            if (model_class == "lm") {
                ps_anova <- prettify(anova(mod))
                nc <- ncol(ps_anova)
                expect_match(ps_anova[, nc - 1], "<0.001")
                expect_match(as.character(ps_anova[, nc]),
                             "\\*\\*\\*")
            }
            if (require("car")) {
                ps_anova <- prettify(Anova(mod))
                ps_anova_lbl <- prettify(Anova(mod),
                                         labels = c(x = "Predictor x",
                                                    z = "nonsense"))
                nc <- ncol(ps_anova)
                if (model_class %in% c("lm", "glm", "lmer", "coxph")) {
                    expect_match(ps_anova[, nc - 1],
                                 "<0.001", info = model_class)
                    expect_match(as.character(ps_anova[, nc]),
                                 "\\*\\*\\*", info = model_class)
                    expect_match(ps_anova_lbl[, 1],
                                 "Predictor x", info = model_class)
                }
                if (model_class == "lme") {
                    expect_match(ps_anova[2, nc - 1],
                                 "<0.001", info = model_class)
                    expect_match(as.character(ps_anova[2, nc]),
                                 "\\*\\*\\*", info = model_class)
                    expect_match(ps_anova_lbl[2, 1],
                                 "Predictor x", info = model_class)
                }
            }
        }
    })
}

################################################################################
## Test HR
################################################################################

if (require("survival")) {
    test_that("survival works", {
        RES <- fit_model("coxph")
        mod <- RES$mod
        data <- RES$data
        CI <- matrix(c(1, 2), ncol = 2)

        ps1 <- prettify(summary(mod), HR = TRUE)
        ps2 <- prettify(summary(mod), HR = FALSE)
        ps3 <- prettify(summary(mod), confint = CI, HR = TRUE)
        ps4 <- prettify(summary(mod), confint = CI, HR = FALSE)

        expect_equivalent(ps1[,2], coef(mod))
        expect_equivalent(ps1[,3], exp(coef(mod)))
        expect_equivalent(ps2[,2], coef(mod))
        expect_equivalent(ps1[,4:5], exp(ps2[, 3:4]))
        expect_equivalent(as.matrix(ps3[,4:5]), exp(CI))
        expect_equivalent(as.matrix(ps4[,3:4]), CI)
    })
}


################################################################################
## Extra column for value labels
################################################################################

fit_anova <- function(model_class =  c("lm", "glm", "coxph", "lme", "lmer")) {
    x <- as.factor(sample(c("Gr. A", "Gr. B", "Gr. C"), 100, replace = TRUE))
    y <- rnorm(100, mean = 2 * (x == "Gr. B") - 1 * (x == "Gr. C"))
    data <- data.frame(y = y, x = x)

    if (model_class %in% c("lme", "lmer")) {
        group <- sample(1:2, 100, replace = TRUE)
        data$group <- group
    }

    if (model_class %in% "coxph") {
        event <- as.logical(sample(0:1, 100, replace = TRUE))
        data$event <- event
        data$y <- exp(y)
    }

    ## fit model with data argument
    mod <- switch(model_class,
                  lm = lm(y ~ x, data = data),
                  glm = glm(y ~ x, data = data),
                  coxph = coxph(Surv(y, event) ~ x, data = data),
                  lme = lme(y ~ x, random = ~ 1 | group, data = data),
                  lmer = lmer(y ~ x + (1 | group), data = data))
    return(list(data = data, model = mod))
}


if (require("nlme") && require("survival") && require("lme4") && packageDescription("lme4")$Version >= 1) {
    test_that("prettify.anova works", {
        for (model_class in c("lm", "glm", "lme", "lmer", "coxph")) {
            ## fit model
            RES <- fit_anova(model_class)
            mod <<- RES$mod
            data <- RES$data

            sum1 <- prettify(summary(mod), confint = FALSE)
            sum2 <- prettify(summary(mod), confint = FALSE, extra.column = TRUE)
            expect_equal(colnames(sum2)[2], "Factor Level", info = model_class)
            if (model_class == "coxph") {
                expect_equal(sum1[1, 1], "x: Gr. B", info = model_class)
                expect_equal(sum2[1, 2], "Gr. B", info = model_class)
                expect_true(all(sum2[,1] == c("x", "")))
            } else {
                expect_equal(sum1[2, 1], "x: Gr. B", info = model_class)
                expect_equal(sum2[2, 2], "Gr. B", info = model_class)
                expect_true(all(sum2[2:3,1] == c("x", "")))
            }
        }
    })
}
