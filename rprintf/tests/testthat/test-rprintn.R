context("rprintn: number-based formatting")

test_that("rprintn works as expected", {
    expect_equivalent(rprintn("{1:s}", "a"), "a")
    expect_equivalent(rprintn("{1:d}", 10), "10")
    expect_equivalent(rprintn("{1:.2f}", pi), "3.14")
    expect_equivalent(rprintn("{1:+.2f}", pi), "+3.14")
    expect_equivalent(rprintn("{1: .2f}", pi), " 3.14")
    expect_equivalent(rprintn("{2:s},{1}", "a", "b"), "b,a")
    expect_equivalent(rprintn("{1},{2:s},{3:.2f},{2},{1}", "a", "b", 1.4), "a,b,1.40,b,a")
    expect_equivalent(rprintn("{{1}}", "a"), "{1}")
    expect_equivalent(rprintn("{{1:s}}", "a"), "{1:s}")
    
    p <- list(name = "Ken", age = 25)
    expect_equivalent(rprintn("{1}, {2}", p), "Ken, 25")
}) 
