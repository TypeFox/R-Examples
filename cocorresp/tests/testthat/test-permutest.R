## Test permutest() method

library("testthat")
library("cocorresp")

context("Testing permutest() method")

## load data
data(beetles, plants, package = "cocorresp")
beetles <- log(beetles + 1)            # log transform the bettle data

test_that("permutest() works", {
    skip_on_cran()
    expect_message(bp.pred <- coca(beetles ~ ., data = plants),
                   regexp = "some species contain no data")
    bp.perm <- permutest(bp.pred, permutations = 499)
    expect_is(bp.perm, "permutest.coca")
    expect_named(bp.perm, c("pval","permstat","total.inertia","inertia",
                            "fitax","pcent.fit","n.axes", "call"))
    expect_output(print(bp.perm),
                  regexp = "Permutation test for predictive co-correspondence analysis:")
})
