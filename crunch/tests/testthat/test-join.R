context("Joining datasets")

join1 <- data.frame(keyvar=c(2, 4, 5, 3), v1=factor(letters[c(2,4)]))
join2 <- data.frame(keyvar=10:1, v2=factor(LETTERS[1:5]))

test_that("join validatation", {
    expect_error(joinDatasets(join1, join2), "x must be a Crunch Dataset")
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(join1, "left"), {
            with(test.dataset(join2, "right"), {
                test_that("Join test setup", {
                    expect_identical(dim(left), dim(join1))
                    expect_identical(dim(right), dim(join2))
                    expect_identical(length(joins(left)), 0L)
                })

                test_that("join validatation", {
                    expect_error(joinDatasets(left, join2),
                        "y must be a Crunch Dataset")
                    expect_error(joinDatasets(left, join2),
                        "y must be a Crunch Dataset")
                    expect_error(joinDatasets(left, right, by=c("v1", "v2")),
                        "Can only join 'by' a single key")
                    expect_error(joinDatasets(left, right, by="v1"),
                        "Join key not found in second Dataset")
                    expect_error(joinDatasets(left, right, by="v2"),
                        "Join key not found in first Dataset")
                })

                joined <- try(joinDatasets(left, right, by="keyvar"))
                test_that("The join succeeded", {
                    expect_true(is.dataset(joined))
                    expect_identical(length(joins(joined)), 1L)
                    skip("TODO: fetch joined variable catalogs")
                    expect_identical(dim(joined), c(4L, 3L))
                    expect_identical(names(joined), c("keyvar", "v1", "v2"))
                    expect_identical(as.vector(joined$v2),
                        factor(c("D", "B", "A", "C")))
                })
            })
        })
    })
}
