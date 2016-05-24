

test_that("Test ScoreDCM outputs for Qmatrix with no interaction", {
  out <- ScoreDCM(observations = observations.test, qmatrix = qmatrix.test, parameter.means = parameter.means.DCM.Mplus.test)
  class_probabilities <- melt(out@results@attribute.profile.result@results, id.vars = "id", variable.name = 'kclass', value.name = 'mcmc_probs')
  test_class_probabilities <- melt(class.probabilities.test, id.vars = "id", variable.name = 'kclass', value.name = 'mplus_probs')
  all_probabilities <- merge(class_probabilities, test_class_probabilities, by = c("id", "kclass"))
  tolerance <- 0.05
  rmse_probs <- sqrt(mean(with(subset(all_probabilities, kclass != 'max.class'), (mcmc_probs - mplus_probs)))^2)
  rmse_class <- sqrt(mean(with(subset(all_probabilities, kclass == 'max.class'), (mcmc_probs - mplus_probs)))^2)
  expect_is(out, "dcm.scorer.class")
  expect_is(out@results, "all.results.class")
  expect_is(out@results@attribute.profile.result, "attribute.profile.class")
  expect_is(out@results@attribute.result, "attribute.class")
  expect_is(out@results@parameter.result, "parameter.class")
  expect_equal(length(out@mcmc.output$class.result), out@mcmc.inputs$nchains)
  expect_less_than(rmse_probs, tolerance)
  expect_less_than(rmse_class, tolerance)    
})



test_that("Test ScoreDCM outputs for Qmatrix with interaction", {
  out <- ScoreDCM(observations = observations.test, qmatrix = qmatrix.test.interaction, parameter.means = parameter.means.DCM.Mplus.interaction.test)
  class_probabilities <- melt(out@results@attribute.profile.result@results, id.vars = "id", variable.name = 'kclass', value.name = 'mcmc_probs')
  test_class_probabilities <- melt(class.probabilities.interaction.test, id.vars = "id", variable.name = 'kclass', value.name = 'mplus_probs')
  all_probabilities <- merge(class_probabilities, test_class_probabilities, by = c("id", "kclass"))
  tolerance <- 0.05
  rmse_probs <- sqrt(mean(with(subset(all_probabilities, kclass != 'max.class'), (mcmc_probs - mplus_probs)))^2)
  rmse_class <- sqrt(mean(with(subset(all_probabilities, kclass == 'max.class'), (mcmc_probs - mplus_probs)))^2)
  expect_is(out, "dcm.scorer.class")
  expect_is(out@results, "all.results.class")
  expect_is(out@results@attribute.profile.result, "attribute.profile.class")
  expect_is(out@results@attribute.result, "attribute.class")
  expect_is(out@results@parameter.result, "parameter.class")
  expect_equal(length(out@mcmc.output$class.result), out@mcmc.inputs$nchains)
  expect_less_than(rmse_probs, tolerance)
  expect_less_than(rmse_class, tolerance)    
})