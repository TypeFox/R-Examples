context("Testing Gaussianize \n")
YY <- matrix(rt(n = 1000, df = 6), ncol = 2)

test_that("has the right attributes", {
  XX <- Gaussianize(YY, type = "h")
  expect_true(sum(grepl("Gaussianized", names(attributes(XX)))) == 4)
})
out <- Gaussianize(YY, type = "h", return.tau.mat = TRUE, 
                   verbose = FALSE, method = "IGMM") 

test_that("inverse works correctly produces identical input", {      
  expect_identical(ncol(out$tau.mat), ncol(YY))
  YY.hat <- Gaussianize(data = out$input, tau.mat = out$tau.mat, 
                        inverse = TRUE)
  expect_equal(lp_norm(YY.hat - YY, 2), 0)
  
  UU <- Gaussianize(YY, type = "h", tau.mat = out$tau.mat, 
                    verbose = FALSE, return.u = TRUE) 
  YY.hat.from.U <- Gaussianize(input.u = UU, tau.mat = out$tau.mat, 
                               inverse = TRUE)
  expect_equal(lp_norm(YY.hat.from.U - YY, 2), 0)
})

test_that("return.u actually returns zero-mean and unit variance input", {
  UU <- Gaussianize(YY, type = "h", tau.mat = out$tau.mat, 
                    verbose = FALSE, return.u = TRUE) 
  expect_equal(lp_norm(colMeans(UU), 1), 0, tol = 1e-4)
  expect_equal(lp_norm(apply(UU, 2, sd) - 1, 1), 0, tol = 1e-3)
})
