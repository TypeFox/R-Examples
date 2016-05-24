context("UNFv6: Numerics")
test_that("Examples from original R package documentation", {
    expect_equal(unf(1:20)$unf, "/FIOZM/29oC3TK/IE52m2A==")
    expect_equal(unf(-3:3)$unf, "7FsSuKWGIp6i7b0NFjckZQ==") # w/ dvn_zero: "pwzm1tdPaqypPWRWDeW6Jw=="
})

test_that("Large values rounded to correct digits", {
    expect_equal(unf6(1111111500)$unf,
                 unf6(1111111600)$unf)
    expect_equal(unf6(1111112400)$unf,
                 unf6(1111112500)$unf)
    expect_equal(unf6(-1111111500)$unf,
                 unf6(-1111111600)$unf)
    expect_equal(unf6(-1111112400)$unf,
                 unf6(-1111112500)$unf)
})

test_that("Small decimal values rounded to correct digits", {
    expect_equal(unf6(.11111115)$unf,
                 unf6(.11111116)$unf)
    expect_equal(unf6(.11111124)$unf,
                 unf6(.11111125)$unf)
    expect_equal(unf6(-.11111115)$unf,
                 unf6(-.11111116)$unf)
    expect_equal(unf6(-.11111124)$unf,
                 unf6(-.11111125)$unf)
})

test_that("Examples from v6 Specification", {
    expect_equal(unf6(c("+1.234568e+",NA,"+0.e+")), unf(c(1.23456789,NA,0)))
    expect_equal(unf(c(1.23456789,NA,0))$unf, "Do5dfAoOOFt4FSj0JcByEw==")
})

test_that("Irrelevant numeric rounding irrelevant", {
    expect_equal(unf6(1)$unf, unf6(1+1e-7)$unf)
    expect_equal(unf6(1, digits = 6)$unf, unf6(1+1e-6, digits=6)$unf)
})

test_that("Relevant numeric rounding relevant", {
    expect_false(unf6(1)$unf == unf6(1+1e-6)$unf)
    expect_false(unf6(1:20)$unf == unf6((1:20) + 1e-5)$unf)
})

test_that("Numeric representations match appropriate character representations", {
    expect_equal(unf6("+1.e+"), unf6(1))
    expect_equal(unf6("+1.e+"), unf6(TRUE))
    expect_equal(unf6("+0.e+"), unf6(0))
    expect_equal(unf6("+0.e+"), unf6(FALSE))
    expect_equal(unf6("-3.e+2"), unf6(-300))
    expect_equal(unf6("+7.3e-4"), unf6(0.00073))
    expect_equal(unf6("+1.234568e+"), unf(1.23456789))
})
