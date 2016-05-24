test_that("TestModularity returns correct results",
          {
            cor.matrix <- RandomMatrix(10)
            rand.hipots <- matrix(sample(c(1, 0), 30, replace=T), 10, 3)
            hip.array <- CreateHipotMatrix(rand.hipots)
            expect_that(sum(laply(hip.array, isSymmetric)), equals(4))
            expect_that(hip.array, is_a("list"))
            expect_that(dim(hip.array[[1]]), equals(c(10, 10)))
            mod.test <- TestModularity(cor.matrix, rand.hipots)
            expect_that(dim(mod.test), equals(c(4, 6)))
            expect_that(mod.test[,"hypothesis"], equals(c("1", "2", "3", "Full Integration")))
            expect_that(colnames(mod.test), equals(c("hypothesis", "Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")))
            expect_true(all(mod.test[1:4, 2:5] >= -1))
            expect_true(all(mod.test[1:4, 2:5] <= 1))
            expect_that(mod.test[,4]/mod.test[,5], equals(mod.test[,6]))
          }
)
test_that("TestModularity returns correct results for Modularity Hypothesis Index",
{
  cov.matrix <- RandomMatrix(10, 1, 1, 10)
  rand.hipots <- matrix(sample(c(1, 0), 30, replace=T), 10, 3)
  hip.array <- CreateHipotMatrix(rand.hipots)
  mod.test <- TestModularity(cov.matrix, rand.hipots, MHI = TRUE)
  expect_that(dim(mod.test), equals(c(4, 6)))
  expect_that(mod.test[,"hypothesis"], equals(c("1", "2", "3", "Full Integration")))
  expect_that(colnames(mod.test), equals(c("hypothesis", "Rsquared", "Probability", "AVG+", "AVG-", "MHI")))
  expect_true(all(mod.test[1:4, 2:5] >= -1))
  expect_true(all(mod.test[1:4, 2:5] <= 1))
  expect_that((mod.test[,4] - mod.test[,5])/CalcICV(cov.matrix), equals(mod.test[,6]))
}
)

