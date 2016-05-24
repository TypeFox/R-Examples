context("Importing Data Tests")



test_that("edgelists of more than 3 columns are not allowed", {
  testEdgelist1 <- data.frame(col1 = letters[1:10], col2 = letters[1:10],
                              col3 = sample(1:10, 10, replace = TRUE),
                              col4 = rnorm(1:10))
  expect_error(as.conflictmat(testEdgelist1),
               "check your raw data: A edgelist should be of either 2 or 3 columns. If it is a win-loss matrix, the column number should be equal to row number.")
})

test_that("factors are not allowed in edgelist", {
  testEdgelist1 <- data.frame(col1 = letters[1:10], col2 = letters[11:20],
                              stringsAsFactors = TRUE)
  expect_warning(as.conflictmat(testEdgelist1))
})

test_that("only a three-column edgelist is allowed if 'weighted = TRUE'", {
  testEdgelist1 <- data.frame(col1 = letters[1:10], col2 = letters[10:1], stringsAsFactors = FALSE)
  expect_error(as.conflictmat(testEdgelist1, weighted = TRUE),
               "Input a matrix or dataframe with three columns, with the third column being Frequency of the interaction")
})

test_that("A warning raised if dyads in weighted edgelist are not unique", {
  testEdgelist2 <- data.frame(col1 = letters[c(1:10, 1)], col2 = letters[c(10:1, 10)],
                              col3 = sample(1:10, 11, replace = TRUE), stringsAsFactors = FALSE)
  expect_warning(as.conflictmat(testEdgelist2, weighted = TRUE),
               "dyads in the weighted edgelist are not unique; the sum of frequencies is taken for duplicated rows.")
})

test_that("the initiator and the recipient should not be the same", {
  testEdgelist1 <- data.frame(col1 = c(letters[1:11], "a"), col2 = c(letters[11:1], "a"), stringsAsFactors = FALSE)
  expect_warning(as.conflictmat(testEdgelist1, weighted = FALSE))
})

test_that("returns conf.mat", {
  testEdgelist1 <- data.frame(col1 = letters[1:10], col2 = letters[10:1], stringsAsFactors = FALSE)
  expect_is(as.conflictmat(testEdgelist1), "conf.mat")
})

test_that("diagonal of the raw win-loss matrix should be zeros", {
  set.seed(1)
  testMatrix1 <- matrix(sample(1:100, 100, TRUE), 10, 10)
  diag(testMatrix1) <- sample(c(0, 1), 10, TRUE, prob = c(0.9, 0.1))
  expect_warning(as.conflictmat(testMatrix1))
})



test_that("outputs are correct", {
   set.seed(1)
   edgelist1 <- data.frame(col1 = sample(letters[1:26], 100, replace = TRUE),
                           col2 = sample(letters[1:26], 100, replace = TRUE),
                          stringsAsFactors = FALSE)
   edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]

#   edgelist1 <- data.frame(x = letters[1:5], y = letters[6:2])

   expect_equal_to_reference(as.conflictmat(edgelist1), file = "asConflictoutput1.rds")
 })


