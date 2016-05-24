context("Combination")

test_that("Combination", {
    I <- iterpc(5, 3, ordered=FALSE)
    getnext(I)
    expect_that(getnext(I)[3], equals(4))
    expect_that(getcurrent(I)[3], equals(4))
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(3, 3, labels=c("a","b","c"), ordered=FALSE)
    getnext(I)
    expect_that(getnext(I)[3], equals(NULL))

    I <- iterpc(3, 2, ordered=FALSE)
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(c(2,2), 2, labels=c("a", "c"), ordered=FALSE)
    getnext(I)
    expect_that(getnext(I)[2], equals("c"))
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(5, 3, ordered=FALSE,replace=TRUE)
    getnext(I)
    expect_that(getnext(I), equals(c(1,1,2)))
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(5, 5, ordered=FALSE)
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(5, 1, ordered=FALSE)
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(5, 5, labels=1:5, ordered=FALSE)
    expect_that(nrow(getall(I)), equals(getlength(I)))

    I <- iterpc(5, 1, labels=1:5, ordered=FALSE)
    expect_that(nrow(getall(I)), equals(getlength(I)))
})
