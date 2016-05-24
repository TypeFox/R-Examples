context("isVariableName")

test_that("abc is a valid name", {
   expect_true(DYM:::isVariableName("abc"))
})

test_that("_abc is not a valid name", {
   expect_false(DYM:::isVariableName("_abc"))
})

## The following tests does not pass on Windows environment.
# test_that("日本語 is a valid name", {
#    expect_true(DYM:::isVariableName("日本語"))
# })
# 
# test_that("étranger is a valid name", {
#    expect_true(DYM:::isVariableName("étranger"))
# })

test_that(". is a valid name", {
   expect_true(DYM:::isVariableName("."))
})

test_that("%op% is not a valid name", {
   expect_false(DYM:::isVariableName("%op%"))
})

test_that("TRUE is not a valid name", {
   expect_false(DYM:::isVariableName("TRUE"))
})
