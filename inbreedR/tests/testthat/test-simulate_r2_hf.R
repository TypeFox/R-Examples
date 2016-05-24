context("simulate_r2_hf")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

test_that("r2 results matrix is computed", {
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, H_nonInb = 0.5, meanF = 0.2, varF = 0.05,
                                       subsets = c(2,4), reps = 10)$estMat), TRUE)
})

test_that("r2 results matrix is computed with snps", {
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, H_nonInb = 0.5, meanF = 0.2, varF = 0.05,
                                       subsets = c(2,4), reps = 10 ,type = "snps")$estMat), TRUE)
})

test_that("subsets have the right range", {
    expect_error(simulate_r2_hf(n_ind = 20, H_nonInb = 0.5, meanF = 0.2, varF = 0.05,
                             subsets = 0, reps = 10))
    expect_error(simulate_r2_hf(n_ind = 20, H_nonInb = 0.5, meanF = 0.2, varF = 0.05,
                             subsets = c(1,2), reps = 10))
    expect_error(simulate_r2_hf(n_ind = 20, H_nonInb = 0.5, meanF = 0.2, varF = 0.05,
                             subsets = NULL, reps = 10))
})