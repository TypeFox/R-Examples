context("Feature Ignore Singletons")

test_that("we can add remove singletons to a models", {
  model <- coal_model(10, 1) + feat_ignore_singletons()
  expect_true(has_ign_singletons(model))
  expect_false(has_ign_singletons(coal_model(10, 1)))
})


test_that("the function to remove singletons works", {
  ssl <- list(create_test_segsites(), create_test_segsites())
  expect_equal(remove_singletons(ssl),
               list(create_test_segsites()[ , -c(1, 3)],
                    create_test_segsites()[ , -c(1, 3)]))

  expect_equal(remove_singletons(list(create_empty_segsites())),
               list(create_empty_segsites()))
})


test_that("singletons are removed from simulations", {
  model <- coal_model(2, 1) +
    feat_mutation(10) +
    feat_ignore_singletons() +
    sumstat_sfs()
  expect_equal(simulate(model)$sfs, 0)
})
