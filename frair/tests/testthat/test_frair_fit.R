library(frair)
data("gammarus")
context("Testing frair_fit()")

test_that("friar_fit understands text and function calls for models", {
    m1 <- frair_fit(eaten~density, data=gammarus, response='rogersII', 
                    start=list(a = 1.2, h = 0.015), fixed=list(T = 1))
    m2 <- frair_fit(eaten~density, data=gammarus, response=rogersII, 
                    start=list(a = 1.2, h = 0.015), fixed=list(T = 1))
    expect_is(m1, 'frfit')
    expect_is(m2, 'frfit')
    expect_equal(coef(m1), coef(m2))
})

test_that("friar_fit errors on unknown models", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='cheese', 
                           start=list(a = 1.2, h = 0.015), fixed=list(T = 1)),
                 "is not a recognised response")
})

test_that("friar_fit errors on complex models", {
    expect_error(frair_fit(eaten~density+dd, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015), fixed=list(T = 1)),
                 "Only simple formulae")
})

test_that("friar_fit errors when start isn't a list", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=c(a = 1.2, h = 0.015), fixed=list(T = 1)),
                 "must be a list containing single, named numeric values")
})

test_that("friar_fit errors when fixed isn't a list", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015), fixed=c(T = 1)),
                 "must be a list containing single, named numeric values")
})

test_that("friar_fit errors on missing start values (a)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(h = 0.015), fixed=list(T = 1)),
                 "The following item is missing: a")
})

test_that("friar_fit errors on missing start values (h)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2), fixed=list(T = 1)),
                 "The following item is missing: h")
})

test_that("friar_fit errors on missing fixed values (T)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015)),
                 "The following item is missing: T")
})

test_that("friar_fit errors when multiple values are missing (h + T)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2)),
                 "The following items are missing: h, T")
})

test_that("friar_fit errors on totally missing start values", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', fixed=list(T = 1)),
                 "You didn't provide starting values")
})

test_that("friar_fit errors on too many individual start values", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = c(1.2, 2), h = 0.015), fixed=list(T = 1)),
                 "must be single, named numeric values")
})

test_that("friar_fit errors on too many individual fixed values", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015), fixed=list(T = c(1,2))),
                 "must be single, named numeric values")
})

test_that("friar_fit errors on an unneccesary start variable (q)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015, q = 2), fixed=list(T = 1)),
                 "The following item is not needed: q")
})

test_that("friar_fit errors on an unneccesary fixed variable (P)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015), fixed=list(T = 1, P = 1)),
                 "The following item is not needed: P")
})

test_that("friar_fit errors on multiple unneccesary start + fixed variables (q + P)", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = 1.2, h = 0.015, q = 2), fixed=list(T = 1, P = 1)),
                 "The following items are not needed: q, P")
})

test_that("friar_fit errors on too many individual start values", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(a = c(1.2, 2), h = 0.015), fixed=list(T = 1)),
                 "must be single, named numeric values")
})

test_that("friar_fit errors on unnamed values", {
    expect_error(frair_fit(eaten~density, data=gammarus, response='rogersII', 
                           start=list(1.2, 0.015), fixed=list(1)),
                 "must be a list containing single, named numeric values")
})