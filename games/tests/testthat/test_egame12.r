library("games")
context("egame12 model")

## Same model as in example
data("war1800")
formula_main <- esc + war ~ s_wt_re1 + revis1 | 0 | regime1 | balanc + regime2
formula_main <- as.Formula(formula_main)
fit_main <- egame12(formula_main, data = war1800)

## Convert to factor outcome
war1800 <- transform(war1800,
                     y = factor(
                         ifelse(esc == 1,
                                ifelse(war == 1,
                                       "war",
                                       "cap"),
                                "sq"),
                         levels = c("sq", "cap", "war")))
formula_fact <- update(formula_main, y ~ .)
fit_fact <- egame12(formula_fact, data = war1800)

test_that("factor- and matrix-form LHS yield same results", {
    expect_equal(coef(fit_main), coef(fit_fact), check.attributes = FALSE)
    expect_equal(vcov(fit_main), vcov(fit_fact), check.attributes = FALSE)
    expect_equal(logLik(fit_main), logLik(fit_fact))
})

test_that("log-likelihood is correctly computed", {
    ## Using `predict` method
    pp <- predict(fit_main)
    pp_observed <- pp[cbind(seq_along(fit_main$y), as.numeric(fit_main$y))]
    expect_equal(fit_main$log.likelihood, log(pp_observed))

    ## Calculated by hand
    cf <- coef(fit_fact)
    u1_sq <- with(fit_fact$model,
                  cf["u1(sq):(Intercept)"]
                  + cf["u1(sq):s_wt_re1"] * s_wt_re1
                  + cf["u1(sq):revis1"] * revis1)
    u1_cap <- 0
    u1_war <- with(fit_fact$model,
                   cf["u1(war):(Intercept)"]
                   + cf["u1(war):regime1dem"] * as.numeric(regime1 == "dem")
                   + cf["u1(war):regime1mixed"] * as.numeric(regime1 == "mixed"))
    u2_war <- with(fit_fact$model,
                   cf["u2(war):(Intercept)"]
                   + cf["u2(war):balanc"] * balanc
                   + cf["u2(war):regime2dem"] * as.numeric(regime2 == "dem")
                   + cf["u2(war):regime2mixed"] * as.numeric(regime2 == "mixed"))
    p4 <- pnorm(u2_war, sd = sqrt(2))
    p3 <- 1 - p4
    p2 <- pnorm(p4 * u1_war + p3 * u1_cap - u1_sq, sd = sqrt(2))
    p1 <- 1 - p2
    p_war <- p2 * p4
    p_cap <- p2 * p3
    p_sq <- p1
    p_observed <- ifelse(fit_fact$y == "sq",
                         p_sq,
                         ifelse(fit_fact$y == "cap",
                                p_cap,
                                p_war))
    expect_equal(log(p_observed), log(pp_observed))
})

test_that("egame12 works with integer-class binary LHS", {
    war1800 <- transform(war1800,
                         esc = as.integer(esc),
                         war = as.integer(war))
    expect_true({
        fit_main_int <- egame12(formula_main, data = war1800)
        TRUE
    })
})
