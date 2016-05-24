context("Feature Sample")

test_that("Creation of sample features works", {
  expect_equal(feat_sample(2)$get_sizes(), 2)
  expect_equal(feat_sample(1:5)$get_sizes(), 1:5)
  expect_error(feat_sample("blub"))
  expect_error(feat_sample(numeric(0)))
})


test_that("sample sizes are reported corrently", {
  expect_equal(get_sample_size(coal_model(c(10, 15))), c(10, 15))
  expect_equal(get_sample_size(coal_model(10)), 10)
})


test_that("sample sizes are reported corrently in polyploid models", {
  expect_equal(get_sample_size(coal_model(c(10, 15), ploidy = 2)), c(20, 30))
  expect_equal(get_sample_size(coal_model(c(10, 15), ploidy = 3)), c(30, 45))
  expect_equal(get_sample_size(coal_model(3, ploidy = 8)), 24)
})


test_that("generating scrm cmd works", {
  model <- coal_model(15, 1)
  expect_equal(get_simulator("scrm")$get_cmd(model), "scrm 15 1 ")
  model <- coal_model(15:16, 1)
  expect_equal(get_simulator("scrm")$get_cmd(model), "scrm 31 1 -I 2 15 16 ")
  model <- coal_model(15:17, 1)
  expect_equal(get_simulator("scrm")$get_cmd(model), "scrm 48 1 -I 3 15 16 17 ")
})
