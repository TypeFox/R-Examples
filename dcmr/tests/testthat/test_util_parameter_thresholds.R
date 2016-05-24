context("Util Parameter Thresholds")

test_that("Test model returns correct threshold parameter values for model DCM", {
  tolerance <- .001
  nattributes <- 3
  nclasses <- 2^nattributes
  pmatrix.test <- GetAttributeProfiles(nattributes)
  threshold.labels <- GetThresholdLabels(qmatrix.test, pmatrix.test)
  parameter.means.DCM.Mplus.test <- c(parameter.means.DCM.Mplus.test, 0)
  names(parameter.means.DCM.Mplus.test) <- c(parameter.means.names.DCM.Mplus.test, "mu_8")
  taus <- parameter.means.DCM.Mplus.test[grep('^tau', names(parameter.means.DCM.Mplus.test))]
  threshold.values <- GetThresholdValues(threshold.labels, taus)
  
  class1 <- c(0.657, 0.138, -0.123, 1.069, 1.423, 0.281, 0.44, 1.955, 0.385, -0.744, 0.842)
  class2 <- c(0.657, 0.138, -2.727, -0.631, 1.423, 0.281, 0.44, 1.955, 0.385, -0.744, 0.842)
  class3 <- c(0.657, -1.247, -0.123, 1.069, 0.396, 0.281, 0.44, 1.955, -1.973, -4.09, -0.551)
  class4 <- c(0.657, -1.247, -2.727, -0.631, 0.396, 0.281, 0.44, 1.955, -1.973, -4.09, -0.551)
  class5 <- c(-0.446, 0.138, -0.123, 1.069, 1.423, -1.522, -2.334, 0.724, 0.385, -0.744, 0.842)
  class6 <- c(-0.446, 0.138, -2.727, -0.631, 1.423, -1.522, -2.334, 0.724, 0.385, -0.744, 0.842)
  class7 <- c(-0.446, -1.247, -0.123, 1.069, 0.396, -1.522, -2.334, 0.724, -1.973, -4.09, -0.551)
  class8 <- c(-0.446, -1.247, -2.727, -0.631, 0.396, -1.522, -2.334, 0.724, -1.973, -4.09, -0.551)
  
  expect_equal(as.vector(threshold.values[1, ]), class1, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[2, ]), class2, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[3, ]), class3, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[4, ]), class4, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[5, ]), class5, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[6, ]), class6, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[7, ]), class7, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[8, ]), class8, tolerance = tolerance)
  
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[1, ])), class1, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[2, ])), class2, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[3, ])), class3, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[4, ])), class4, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[5, ])), class5, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[6, ])), class6, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[7, ])), class7, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus, threshold.labels[8, ])), class8, tolerance = tolerance)
  
  # DCM model with three attrs and one interaction in qmatrix
  threshold.labels <- GetThresholdLabels(qmatrix.test.interaction, pmatrix.test)
  parameter.means.DCM.Mplus.interaction.test <- c(parameter.means.DCM.Mplus.interaction.test, 0)
  names(parameter.means.DCM.Mplus.interaction.test) <- c(parameter.means.names.DCM.Mplus.interaction.test, "mu_8")
  taus.interaction <- parameter.means.DCM.Mplus.interaction.test[grep('^tau', names(parameter.means.DCM.Mplus.interaction.test))]
  threshold.values <- GetThresholdValues(threshold.labels, taus.interaction)

  class1.interaction = c(0.648, 0.065, 0.015, 1.111, 1.415, 0.235, 0.285, 2.011, 0.61, -0.647, 0.952)
  class2.interaction = c(0.648, 0.065, -2.695, -0.611, 1.415, 0.235, 0.285, 2.011, 0.61, -0.647, 0.952)
  class3.interaction = c(0.648, 0.0649, 0.015, 1.111, 0.413, 0.235, 0.285, 2.011, -1.987, -3.996, -0.542)
  class4.interaction = c(0.648, 0.0649, -2.695, -0.611, 0.413, 0.235, 0.285, 2.011, -1.987, -3.996, -0.542)
  class5.interaction = c(-0.482, -1.414, 0.015, 1.111, 1.415, -1.587, -2.39, 0.684, 0.61, -0.647, 0.952)
  class6.interaction = c(-0.482, -1.414, -2.695, -0.611, 1.415, -1.587, -2.39, 0.684, 0.61, -0.647, 0.952)
  class7.interaction = c(-0.482, -1.4241, 0.015, 1.111, 0.413, -1.587, -2.39, 0.684, -1.987, -3.996, -0.542)
  class8.interaction = c(-0.482, -1.4241, -2.695, -0.611, 0.413, -1.587, -2.39, 0.684, -1.987, -3.996, -0.542)
  
  expect_equal(as.vector(threshold.values[1, ]), class1.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[2, ]), class2.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[3, ]), class3.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[4, ]), class4.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[5, ]), class5.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[6, ]), class6.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[7, ]), class7.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[8, ]), class8.interaction, tolerance = tolerance)
  
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[1, ])), class1.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[2, ])), class2.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[3, ])), class3.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[4, ])), class4.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[5, ])), class5.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[6, ])), class6.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[7, ])), class7.interaction, tolerance = tolerance)
  expect_equal(as.vector(GetClassThresholds(taus.interaction, threshold.labels[8, ])), class8.interaction, tolerance = tolerance)
  
  # Testing for Kernel values 
  # DCM
  model.type <- "DCM"
  threshold.labels <- GetThresholdLabels(qmatrix.test.interaction, pmatrix.test)
  lambda.equations <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)$lambdas$lambda.equations
  names(parameter.means.DCM.kernel.Mplus.interaction.test) <- parameter.means.names.DCM.kernel.Mplus.interaction.test
  threshold.info <- GetThresholdValuesKernel(lambda.equations, parameter.means.DCM.kernel.Mplus.interaction.test, threshold.labels)
  taus.interaction.from.kernel <- -threshold.info$taus
  threshold.values <- -threshold.info$threshold.values
  
  expect_equal(taus.interaction.from.kernel, taus.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[1, ]), class1.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[2, ]), class2.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[3, ]), class3.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[4, ]), class4.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[5, ]), class5.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[6, ]), class6.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[7, ]), class7.interaction, tolerance = tolerance)
  expect_equal(as.vector(threshold.values[8, ]), class8.interaction, tolerance = tolerance)

  # NCRUM
  model.type <- "NCRUM"
  threshold.labels <- GetThresholdLabels(qmatrix.test.interaction, pmatrix.test)
  # Interaction (2 rias for item with interaction)
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  taus.ncrum.population.values <- c(-1.798, -2.991, -2.221, -1.86, -2.599, -1.067, -3.009, -0.852, -1.858
                                    , -2.643, -1.044, 1.802, 2.865, -0.666, 1.086, 0.706, 0.442, 3.013, 1.465
                                    , -1.92, 2.893, 0.701, 2.377, 2.284)
  
  names(parameter.means.NCRUM.interaction.test) <- parameter.means.names.NCRUM.interaction.test
  taus.ncrum <- GetThresholdValuesKernelPiR(parameter.means.NCRUM.interaction.test, qmatrix.test.interaction, pmatrix.test, threshold.labels)
  expect_equal(as.vector(taus.ncrum), taus.ncrum.population.values, tolerance = tolerance)
})