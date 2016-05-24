
context("duplicated")
test_that("duplicated.editmatrix works",{

    expect_equal(
        duplicated.editmatrix(editmatrix(expression(
            x < 0, x <= 0, y < 0))),
        c(FALSE,FALSE,FALSE)
    )
})



