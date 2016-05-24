context("check_if_data_matrix")

test_that("check_if_data_matrix behaves as expected", {
    expect_true(check_if_data_matrix(matrix(0)))
    expect_true(check_if_data_matrix(data.frame()))
    expect_false(check_if_data_matrix(list()))
    expect_false(check_if_data_matrix(numeric(10)))

    ### Check Matrix package
# 4-12-16: Why do we want Matrix objects to return FALSE here? Deprecated for now.
#     library("Matrix")
#     expect_false(check_if_data_matrix(Matrix::Diagonal(10, 1)))
#     expect_false(check_if_data_matrix(Matrix::Matrix(matrix(0))))

})
