context("reIndex")

test_that("reIndexC <-> reIndexR returns original", {
    ### R -> C
    sp <- sparse(list(rows = c(1L), cols = c(1L), vals = pi, dim = c(2, 2), start = 1))
    expect_equal(reIndexR(reIndexC(sp)), sp)

#     sbm <- SparseBlockMatrixR.list(list(rows = list(c(1L)), vals = list(c(pi)), blocks = list(integer(0)), sigmas = exp(1), start = 1))
#     expect_equal(reIndexR(reIndexC(sbm)), sbm)

    ### C -> R
    sp[["start"]] <- 0 # pretend like this is uses C-style indexing
    expect_equal(reIndexC(reIndexR(sp)), sp)

#     sbm[["start"]] <- 0 # pretend like this is uses C-style indexing
#     expect_equal(reIndexC(reIndexR(sbm)), sbm)
})

test_that("reIndexC works as expected", {
    ### test on sparse
    spR <- sparse(list(rows = c(1L), cols = c(1L), vals = pi, dim = c(2, 2), start = 1))
    spC <- reIndexC(spR)
    expect_equal(spC, sparse.list(list(rows = c(0L), cols = c(0L), vals = pi, dim = c(2, 2), start = 0)))

    ### test on SparseBlockMatrixR
#     sbmR <- SparseBlockMatrixR.list(list(rows = list(c(1L)), vals = list(c(pi)), blocks = list(integer(0)), sigmas = exp(1), start = 1))
#     sbmC <- reIndexC(sbmR)
#     expect_equal(sbmC, SparseBlockMatrixR.list(list(rows = list(c(0L)), vals = list(c(pi)), blocks = list(integer(0)), sigmas = exp(1), start = 0)))
})

test_that("reIndexR works as expected", {
    ### test on sparse
    spC <- sparse(list(rows = c(0L), cols = c(0L), vals = pi, dim = c(2, 2), start = 0))
    spR <- reIndexR(spC)
    expect_equal(spR, sparse.list(list(rows = c(1L), cols = c(1L), vals = pi, dim = c(2, 2), start = 1)))

    ### test on SparseBlockMatrixR
#     sbmC <- SparseBlockMatrixR.list(list(rows = list(c(0L)), vals = list(c(pi)), blocks = list(integer(0)), sigmas = exp(1), start = 0))
#     sbmR <- reIndexR(sbmC)
#     expect_equal(sbmR, SparseBlockMatrixR.list(list(rows = list(c(1L)), vals = list(c(pi)), blocks = list(integer(0)), sigmas = exp(1), start = 1)))
})
