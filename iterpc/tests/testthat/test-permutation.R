context("Permutation")

test_that("Permutation", {
    I <- iterpc(5, ordered=TRUE)
    expect_that(getnext(I)[3], equals(3))
    expect_that(getcurrent(I)[3], equals(3))
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(3, labels=c("a","b","c"), ordered=TRUE)
    getnext(I)
    expect_that(getnext(I), equals(c("a","c","b")))
    expect_that(getcurrent(I), equals(c("a","c","b")))

    I <- iterpc(c(2,1), labels=c("a","c"), ordered=TRUE)
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(c(2,2), 2, labels=c("a","c"), ordered=TRUE)
    expect_that(nrow(getall(I)), equals(getlength(I)))
    getnext(I)
    expect_that(getnext(I)[2], equals("c"))

    I <- iterpc(c(2,2), 2, labels=c("a","c"), replace=TRUE, ordered=TRUE)
    expect_that(nrow(getall(I)), equals(getlength(I)))
    expect_that(getnext(I,2)[2,2], equals("c"))

    I <- iterpc(5, 1, ordered=TRUE)
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(5, 1, labels=1:5, ordered=TRUE)
    expect_that(nrow(getall(I)), equals(getlength(I)))
})
