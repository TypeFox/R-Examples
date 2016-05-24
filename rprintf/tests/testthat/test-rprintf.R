context("rprintf: general formatting")

test_that("rprintf works as expected", {
    expect_equivalent(rprintf("%s", "a"), "a")
    expect_equivalent(rprintf("$name:s", name = "a"), "a")
    expect_equivalent(rprintf("{1:s}", "a"), "a")
    expect_equivalent(rprintf("$name", name = "a"), "a")
    expect_equivalent(rprintf("{1}", "a"), "a")
    
    p <- list(name = "Ken", age = 25)
    expect_equivalent(rprintf("$name, $age", p), "Ken, 25")
}) 
