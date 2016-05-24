test_that("DriftTest returns resonable results",
{
  means <- array(rnorm(40*10), c(10, 40)) 
  means.list <- alply(means, 1)
  cov.matrix <- RandomMatrix(40, 1, 1, 10)
  test.array <- DriftTest(means, cov.matrix, FALSE)
  test.list <- DriftTest(means.list, cov.matrix, FALSE)
  expect_equal(test.array[-6], test.list[-6])
  expect_is(test.array, "list")
  expect_is(test.array[[1]], "lm")
  expect_is(test.array[[2]], "matrix")
  expect_is(test.array[[3]], "numeric")
  expect_is(test.array[[4]], "numeric")
  expect_is(test.array[[5]], "logical")
  expect_is(test.array[[6]], c("gg", "ggplot"))
  expect_equal(names(test.array), c("regression",
                                    "coefficient_CI_95",
                                    "log.between_group_variance",
                                    "log.W_eVals",
                                    "drift_rejected",
                                    "plot"))
  expect_equal(length(test.array$log.between_group_variance), dim(means)[2])
  expect_equal(length(test.array$log.W_eVals), dim(means)[2])
  expect_equal(test.array$log.W_eVals, log(eigen(cov.matrix)$values))
})

test_that("TreeDriftTest returns resonable results",
{
  library(ape)
  data(bird.orders)
  tree <- bird.orders
  mean.list <- llply(tree$tip.label, function(x) rnorm(5))
  names(mean.list) <- tree$tip.label
  cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
  names(cov.matrix.list) <- tree$tip.label
  w.cov <- PhyloW(tree, cov.matrix.list)$'24'
  test.list <- TreeDriftTest(tree, mean.list, cov.matrix.list)
  expect_is(test.list, 'list')
  expect_equal(length(test.list), 13)
  expect_equal(test.list[[length(test.list)]][-6], 
               DriftTest(mean.list, w.cov, FALSE)[-6])
})

test_that("TreeDriftTest returns resonable results",
          {
            library(ape)
            data(bird.orders)
            tree <- bird.orders
            mean.list <- llply(tree$tip.label, function(x) rnorm(5))
            names(mean.list) <- tree$tip.label
            cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
            names(cov.matrix.list) <- tree$tip.label
            cov.matrix.list = cov.matrix.list[sample(1:length(cov.matrix.list))]
            mean.list_2 = mean.list[sample(1:length(mean.list))]
            w.cov <- PhyloW(tree, cov.matrix.list)$'32'
            test.list <- TreeDriftTest(tree, mean.list_2, cov.matrix.list)
            expect_is(test.list, 'list')
            expect_equal(length(test.list), 13)
            expect_equal(test.list[["32"]][-6], 
                         DriftTest(mean.list[na.omit(tree$tip.label[phytools::getDescendants(tree, "32")])], 
                                   w.cov, FALSE)[-6])
          })

          