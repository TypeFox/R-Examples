context("Test disparity filter")

test_that("undirected network", {
  
  g <- sample_pa(n = 250, m = 5, directed = FALSE)
  E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)
  backbone(g)

})

test_that("directed network", {

  g <- sample_pa(n = 250, m = 5, directed = TRUE)
  E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)
  backbone(g)
  
})

test_that("errors", {

  g <- data.frame(i = 1:3, j = 2:4, w = sample(1:10, 3))
  expect_error(backbone(g), "Not a graph object")
  
})
