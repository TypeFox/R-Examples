context("Hold-out")

test_that("Hold-out method",{
    V1 <- 1:10
    V2 <- 11:20
    a <- data.frame(V1,V2)
    expect_that(is.function(b <- hold_out(a)), is_true())
    o1 <- data.frame(V1=1,V2=11)
    r1 <- data.frame(V1=2:10,V2=12:20,row.names=2:10)
    res1 <- list(one=o1, rest=r1)
    expect_that(b(), equals(res1))
    for(i in 2:9) b()
    o2 <- a[10,]
    r2 <- data.frame(V1=1:9,V2=11:19,row.names=1:9)
    res2 <- list(one=o2, rest=r2)
    expect_that(b(), equals(res2))
    expect_that(is.na(b()), is_true())
    expect_that(is.na(b()), is_true())
    expect_that(is.na(b()), is_true())
})
