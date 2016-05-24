context("Testing for Starting Points")

test_that("Testing Starting Points - 3D simplex", {
  
  set.seed(314)
  
  ## Initialize 3D Simplex
  
  n <- 3
  A <- matrix(rep(1,n), ncol = n, nrow = 1)
  b <- 1
  
  ## find the particular and homogeneous solutions

  z <- walkr:::complete_solution(A,b)
  particular  <- z$particular
  homogeneous <- z$homogeneous
  
  ## now, the non-negativity constraints could be expressed as
  ## particular + homogeneous %*% alpha >= 0
  ## -homogeneous %*% alpha <= particular
  ## thus, new_A and new_b are as follows
  
  new_A <- -homogeneous
  new_b <- particular
  
  ## note that these sampled points are alphas 
  
  start.points1 <- walkr:::start_point(A = new_A, b = new_b, n = 10, average = 10)
  start.points2 <- walkr:::start_point(A = new_A, b = new_b, n = 5, average = 20)
  
  ## now we use a mapping function to obtain the points on the simplex 
  
  answer1 <- apply(start.points1, 2, function(x) { homogeneous %*% x + particular  })
  answer2 <- apply(start.points2, 2, function(x) { homogeneous %*% x + particular  })
  
  ## no points with elements < 0
  
  expect_equal(length(which(answer1 <0)), 0)
  expect_equal(length(which(answer2 <0)), 0)
  
  ## all points sum up to one
  
  sum1 <- apply(answer1, 2, sum)
  sum2 <- apply(answer2, 2, sum)
  
  expect_true(all(sum1 - 1 <= 1e-10))
  expect_true(all(sum2 - 1 <= 1e-10))
 
})

test_that("Testing Starting Points - 20D with constraints", {
  
  set.seed(314)
  
  ## Initialize 20D Simplex
  
  n <- 20
  A <- matrix(rep(1,n), ncol = n, nrow = 1)
  A <- rbind(A, sample(c(1,0), n, replace = T))
  A <- rbind(A, sample(c(1,0), n, replace = T))
  b <- c(1, 0.5, 0.1)
  
  ## find the particular and homogeneous solutions
  
  z <- walkr:::complete_solution(A,b)
  particular  <- z$particular
  homogeneous <- z$homogeneous
  
  ## now, the non-negativity constraints could be expressed as
  ## particular + homogeneous %*% alpha >= 0
  ## -homogeneous %*% alpha <= particular
  ## thus, new_A and new_b are as follows
  
  new_A <- -homogeneous
  new_b <- particular
  
  ## note that these sampled points are alphas 
  
  start.points1 <- walkr:::start_point(A = new_A, b = new_b, n = 10, average = 10)
  start.points2 <- walkr:::start_point(A = new_A, b = new_b, n = 5, average = 20)
  
  ## now we use a mapping function to obtain the points on the simplex 
  
  answer1 <- apply(start.points1, 2, function(x) { homogeneous %*% x + particular  })
  answer2 <- apply(start.points2, 2, function(x) { homogeneous %*% x + particular  })
  
  ## SIMPLEX: no points with elements < 0
  
  expect_equal(length(which(answer1 <0)), 0)
  expect_equal(length(which(answer2 <0)), 0)
  
  ## SIMPLEX: all points sum up to one
  
  sum1 <- apply(answer1, 2, sum)
  sum2 <- apply(answer2, 2, sum)
  
  expect_true(all(sum1 - 1 <= 1e-10))
  expect_true(all(sum2 - 1 <= 1e-10))
  
  ## CONSTRAINT: all satisfied 
  
  expect_true(all(apply(answer1, 2, function(x) { A %*% x - b < 1e-10  })))
  
  
})
