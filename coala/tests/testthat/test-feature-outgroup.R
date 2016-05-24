context("Feature Outgroup")

test_that("Creation of outgroup feature works", {
  expect_equal(feat_outgroup(1)$get_population(), 1)
  expect_equal(feat_outgroup(2)$get_population(), 2)
})


test_that("Outgroup setting and getting works", {
  model <- coal_model(1:4 * 2, 100)
  expect_equal(get_outgroup(model), NA)

  for (i in 1:4) {
    expect_equal(get_outgroup(model + feat_outgroup(i)), i)
    expect_equal(get_outgroup_size(model + feat_outgroup(i)), 2 * i)
  }
})
