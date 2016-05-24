VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE, colClasses = 'numeric'))
SEQUENCE <- seq_len(37) - 1

context("Testing ddelap")
test_that("Singleton ddelap functions", {
  expect_equal(ddelap(SEQUENCE, 1, 4, 2), VAL$DDELAP_1)
  expect_equal(ddelap(SEQUENCE, 5, 3, 8, log = TRUE), log(ddelap(SEQUENCE, 5, 3, 8)))
})
test_that("Vector ddelap functions", {
  expect_equal(ddelap(SEQUENCE, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)), VAL$DDELAP_Triple)
  expect_equal(ddelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), log = TRUE), 
              log(ddelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))))
})
  
context("Testing pdelap")
test_that("Singleton pdelap functions", {
  expect_equal(pdelap(SEQUENCE, 2, 1, 5), VAL$PDELAP_2)
  expect_equal(pdelap(SEQUENCE, 4, 5, 1, log.p = TRUE), log(pdelap(SEQUENCE, 4, 5, 1)))
  expect_equal(pdelap(seq_len(101) - 1, 8, 10, 6, lower.tail = FALSE), 1 - pdelap(seq_len(101) - 1, 8, 10, 6))
  expect_equal(pdelap(6, 2.9647, 0.005/2.9647, 0.0057, lower.tail=FALSE), 0)
})

test_that("Vector pdelap functions", {
  expect_equal(pdelap(SEQUENCE, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)), VAL$PDELAP_Triple)
  expect_equal(pdelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), lower.tail = FALSE),
              1 - pdelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9)))
})

context("Testing qdelap")
test_that("Singleton qdelap functions", {
  expect_equal(qdelap(.4, 1, 4, 2), 4)
  expect_equal(qdelap(0.4971031395336245, 4, 6, 3, lower.tail = TRUE), 25)
  expect_equal(qdelap(-0.255, 20, 15, 50, log.p = TRUE), 400)
})
test_that("Vector qdelap functions", {
  expect_equal(qdelap(c(.4, .07), c(1, 2), c(4, 1), c(2, 5)), c(4, 3))
  expect_equal(qdelap(c(0.4971031395336245, 0.1374992163369109), c(4, 1), c(6, 9), c(3, 12), lower.tail = FALSE), c(25, 31))
})

context("Testing rdelap")
test_that("Singleton rdelap functions", {
  DP1 <- rdelap(1e6, alpha = 10, beta = 2, lambda = 10) 
  DP2 <- rdelap(1e6, alpha = 2, beta = 14, lambda = 2, exact = FALSE)
  expect_that(abs((mean(DP1) / 30 - 1)) < 5e-3, is_true())
  expect_that(abs((mean(DP2) / 30 - 1)) < 5e-3, is_true())
})
test_that("Vector rdelap functions", {
  DP3 <- rdelap(1e3, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5))
  DP4 <- rdelap(1e6, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5), exact = FALSE)
  expect_that(abs((mean(DP3) / 30 - 1)) < 5e-2, is_true())
  expect_that(abs((mean(DP4) / 30 - 1)) < 5e-3, is_true())
})

context("Testing MoMdelap")
test_that("MoMdelap functions", {
  TestData <- c(5, 7, 9, 9, 10, 11, 11, 13, 17, 24)
  MoM <- MoMdelap(TestData)
  expect_equal(MoM[[1]], 0.769626413810781)
  expect_equal(MoM[[2]], 4.83611111111112)
  expect_equal(MoM[[3]], 7.87800114876508)
})
test_that("MoMdelap traps bad parameters", {
  TestData <- c(3,  2, 12, 11,  1,  7,  1,  4,  0, 4)
  expect_that(MoMdelap(TestData), throws_error("Method of moments not appropriate for this data; results include negative parameters."))
})