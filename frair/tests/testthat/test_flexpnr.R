library(frair)
data("gammarus")
data("bythotrephes")
context("Testing flexpnr")

data(gammarus)
pulex <- gammarus[gammarus$spp=='G.pulex', ]

test_that("friar_fit understands a flexpnr response", {
    m1 <- frair_fit(eaten~density, data=pulex, 
                    response='flexpnr', start=list(b = 1.2, q = 0, h = 0.015), 
                    fixed=list(T=1))
    expect_is(m1, 'frfit')
})

test_that("flexpnr collapses to rogerII when q = 0", {
    expfit <- frair_fit(eaten~density, data = pulex, 
                        response = 'flexpnr', start = list(b = 1.2, h = 0.015), 
                        fixed = list(T = 1, q = 0))
    rogfit <- frair_fit(eaten~density, data=pulex, 
                        response = 'rogersII', start = list(a = 1.2, h = 0.015), 
                        fixed=list(T = 1))
    # NB: Use expect_equivalent for check.attributes = FALSE
    expect_equivalent(coef(expfit)['b'], coef(rogfit)['a']) 
    expect_equal(coef(expfit)['h'], coef(rogfit)['h'])
})


