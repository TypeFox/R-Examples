context("CJ.dt")
library("data.table")
test_that("CJ.dt works correctly on two data.tables", {
    dt.a <- data.table(b = 1:2, c = letters[1:2])
    dt.b <- data.table(d = 3:4, e = letters[3:4])
    ab <- CJ.dt(dt.a, dt.b)
    expect_true(nrow(ab) == nrow(dt.a) * nrow(dt.b))
})

test_that("CJ.dt errors correctly on a data.frame", {
    dt.a <- data.frame(b = 1:2, c = letters[1:2])
    dt.b <- data.table(d = 3:4, e = letters[3:4])
    expect_error(CJ.dt(dt.a, dt.b))
}) 
