test_that("Frechet ES.",{
  # Success - 1
  expect_equal(as.matrix(623.485), FrechetES(3.5, 2.3, 1.6, 10, .95, 30), tolerance=0.01)
  
  # Success - 2
  expect_equal(as.matrix(0.1514), FrechetES(-1.3, 1.2, .2, 10, .9, 280), tolerance=0.01)
  
})