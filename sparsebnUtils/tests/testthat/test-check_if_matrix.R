context("check_if_matrix")

test_that("check_if_matrix behaves as expected", {
    expect_true(check_if_matrix(matrix(0)))
    expect_false(check_if_matrix(list()))
    expect_false(check_if_matrix(numeric(10)))

    ### Check Matrix package
    library("Matrix")
    expect_true(check_if_matrix(Matrix::Diagonal(10, 1)))
    expect_true(check_if_matrix(Matrix::Matrix(matrix(0))))

})
