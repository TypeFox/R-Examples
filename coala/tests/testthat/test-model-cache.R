context("Model Cache")

test_that("get id works", {
  id <- as.numeric(get_id())
  expect_equal(get_id(), as.character(id + 1))
  expect_equal(get_id(), as.character(id + 2))
  expect_equal(get_id(), as.character(id + 3))
})


test_that("caching values works", {
  model1 <- coal_model(5:6, 100)
  cache(model1, "a", 1)
  expect_equal(read_cache(model1, "a"), 1)
  cache(model1, "b", 2)
  expect_equal(read_cache(model1, "b"), 2)

  model2 <- coal_model(5:6, 100)
  cache(model2, "a", 2)
  expect_equal(read_cache(model1, "a"), 1)
  expect_equal(read_cache(model2, "a"), 2)
  expect_equal(read_cache(model2, "b"), NULL)

  model3 <- model2 + feat_mutation(par_const(5))
  cache(model3, "a", 3)
  expect_equal(read_cache(model1, "a"), 1)
  expect_equal(read_cache(model2, "a"), 2)
  expect_equal(read_cache(model3, "a"), 3)
})


test_that("resetting the cache works", {
  model1 <- coal_model(5:6, 100)
  cache(model1, "a", 1)
  expect_equal(read_cache(model1, "a"), 1)
  reset_cache()
  expect_equal(read_cache(model1, "a"), NULL)

  # Check that cache is still usable afterwards
  cache(model1, "a", 2)
  expect_equal(read_cache(model1, "a"), 2)

  model2 <- coal_model(5:6, 100)
  cache(model2, "a", 2)
  expect_equal(read_cache(model2, "a"), 2)
})

