context("Testing quadratic_format()\n")

nn <- 10
kVec <- runif(nn)
kAllZeroVec <- rep(0, nn)

# quadratic matrix of with appropriate dimension to kVec
kMat <- matrix(rnorm(nn^2), ncol = nn)
# this matrix is positive definite
kPosDefMat <- t(Conj(kMat)) %*% kMat

eigen.results <- eigen(kPosDefMat, symmetric = TRUE)
eigen.nonpos <- eigen(kMat)

test_that("quadratic form actually computes the quadratic form", {
  expect_equal(quadratic_form(kMat, kVec), 
               c(Conj(t(kVec)) %*% kMat %*% kVec))
  # 0 times 0 is zero
  expect_equal(0, quadratic_form(matrix(0), 0))
  # first argument must be a matrix
  expect_error(quadratic_form(kVec, kVec))
  # for positive def matrix the quadratic form is positive
  expect_true(quadratic_form(kPosDefMat, kVec) > 0)
  # for eigenvalue it must hold that x' A x = lambda x'x = lambda for
  # unit norm x
  expect_equal(quadratic_form(kPosDefMat, eigen.results$vector[, 1]), 
               eigen.results$value[1])
  
  # for complex values as well
  qf.complex <- quadratic_form(kMat, eigen.nonpos$vector[, 1]) + 0i
  expect_equal(Re(qf.complex),  
               Re(eigen.nonpos$value[1]),
               info = "real part is not equal")
  
  expect_equal(abs(Im(qf.complex)),  
               abs(Im(eigen.nonpos$value[1])),
               info = "imaginary part does not have same magnitude (ignoring sign)")
  
})

context("Testing fill_hermitian()\n")

kRealMat <- matrix(seq_len(16), ncol = 4)
kComplexMat <- kRealMat + 1i * 2 * kRealMat

test_that("fill_hermitian only allows matrices with a real-valued diagonal", {
  expect_error(fill_hermitian(matrix(1 + 1i)))
  expect_error(fill_hermitian(kComplexMat))
})

diag(kComplexMat) <- seq_len(ncol(kComplexMat))

test_that("fill_hermitian only allows NA in lower triangual matrices", {
  expect_error(fill_hermitian(kRealMat))
  expect_error(fill_hermitian(kComplexMat))
  
  kRealMat[lower.tri(kRealMat)] <- 0
  expect_error(fill_hermitian(kRealMat))
  
  kComplexMat[lower.tri(kComplexMat)] <- NA
  expect_true(inherits(fill_hermitian(kComplexMat), "matrix"))
})

kRealMat[lower.tri(kRealMat)] <- NA
kComplexMat[lower.tri(kComplexMat)] <- NA

test_that("fill_hermitian actually produces hermitian matrix", {
  expect_error(fill_hermitian(1))
  expect_equal(fill_hermitian(matrix(1)), matrix(1))
  
  filled.real <- fill_hermitian(kRealMat)
  filled.complex <- fill_hermitian(kComplexMat)
  
  # upper triangual stays the same
  expect_equal(filled.real[upper.tri(filled.real)],
               kRealMat[upper.tri(kRealMat)])
  expect_equal(filled.complex[upper.tri(filled.complex)],
               kComplexMat[upper.tri(kComplexMat)])

  # diagonal stays the same (not doubles)
  expect_equal(diag(filled.real), diag(kRealMat))
  expect_equal(diag(filled.complex), diag(kComplexMat))
  
  # it is actually Hermitian: A = Conj(A)'
  expect_equal(filled.real, t(Conj(filled.real)))
  expect_equal(filled.complex, t(Conj(filled.complex)))
  
})
