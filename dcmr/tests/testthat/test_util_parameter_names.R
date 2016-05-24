context("Util Parameter Names")

test_that("Test model returns correct parameter names (no kernel)", {
  nattributes <- 3
  nclasses <- 2^nattributes
  pmatrix.test <- matrix(as.numeric(unlist(strsplit(binary(0:(nclasses-1)), "")))
                         , nclasses, nattributes, byrow = T)
  parameterization.method <- ''
  # DCM model with three attrs and no interactions in qmatrix
  test.names <- GetParameterNames(qmatrix.test, nattributes, parameterization.method)
  expect_equal(test.names, c(parameter.means.names.DCM.Mplus.test, "mu_8"))
  
  # DCM model with three attrs and one interaction in qmatrix
  test.names <- GetParameterNames(qmatrix.test.interaction, nattributes, parameterization.method)
  expect_equal(test.names, c(parameter.means.names.DCM.Mplus.interaction.test, "mu_8"))
  
  # DCM model with two attrs and no interaction in qmatrix
  nattributes <- 2
  nclasses <- 2^nattributes
  pmatrix.test.attr2 <- matrix(as.numeric(unlist(strsplit(binary(0:(nclasses-1)), "")))
                               , nclasses, nattributes, byrow = T)
  test.names <- GetParameterNames(qmatrix.test, nattributes, parameterization.method)
  parameter.means.names.attr2 <- c("tau1_0", "tau2_0", "tau3_0", "tau4_0", "tau5_0", "tau6_0", "tau7_0"
                                   , "tau8_0", "tau9_0", "tau10_0", "tau11_0", "tau2_2", "tau5_2", "tau9_2"
                                   , "tau10_2", "tau11_2", "tau1_1", "tau6_1", "tau7_1", "tau8_1", "mu_1"
                                   , "mu_2", "mu_3", "mu_4")
  expect_equal(test.names, parameter.means.names.attr2)
  
  nattributes <- 3
  # with parameterization.method = 'Mplus'
  parameterization.method <- 'Mplus'
  # DCM model with three attrs and no interactions in qmatrix
  test.names <- GetParameterNames(qmatrix.test, nattributes, parameterization.method)
  expect_equal(test.names, parameter.means.names.DCM.Mplus.test)
  
  # DCM model with three attrs and one interaction in qmatrix
  test.names <- GetParameterNames(qmatrix.test.interaction, nattributes, parameterization.method)
  expect_equal(test.names, parameter.means.names.DCM.Mplus.interaction.test)
})

test_that("Test model with kernel parameters returns correct parameter names", {
  nattributes <- 3
  nclasses <- 2^nattributes
  pmatrix.test <- matrix(as.numeric(unlist(strsplit(binary(0:(nclasses-1)), "")))
                         , nclasses, nattributes, byrow = T)
  
  # DCM
  model.type <- "DCM"
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  gammas.attr3 <- c('g_1_1', 'g_1_2', 'g_1_3', 'g_2_12', 'g_2_13', 'g_2_23', 'g_3_123')
  lambdas <- c('l1_0', 'l1_1_1', 'l2_0', 'l2_1_2', 'l3_0'
               , 'l3_1_3', 'l4_0', 'l4_1_3', 'l5_0', 'l5_1_2', 'l6_0', 'l6_1_1'
               , 'l7_0', 'l7_1_1', 'l8_0', 'l8_1_1', 'l9_0', 'l9_1_2', 'l10_0'
               , 'l10_1_2', 'l11_0', 'l11_1_2')
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas)
  
  # Testing with an interaction
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  lambdas.interaction <- c('l1_0', 'l1_1_1', 'l2_0', 'l2_1_2', 'l2_1_1', 'l2_2_12', 'l3_0'
               , 'l3_1_3', 'l4_0', 'l4_1_3', 'l5_0', 'l5_1_2', 'l6_0', 'l6_1_1'
               , 'l7_0', 'l7_1_1', 'l8_0', 'l8_1_1', 'l9_0', 'l9_1_2', 'l10_0'
               , 'l10_1_2', 'l11_0', 'l11_1_2')
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.interaction)
  
  # DCM model with two attrs and no interaction in qmatrix
  nattributes <- 2
  nclasses <- 2^nattributes
  pmatrix.test.attr2 <- matrix(as.numeric(unlist(strsplit(binary(0:(nclasses-1)), "")))
                               , nclasses, nattributes, byrow = T)
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  gammas.attr2 <- c('g_1_1', 'g_1_2', 'g_2_12')
  lambdas.attr2.DCM <- c('l1_0', 'l1_1_1', 'l2_0', 'l2_1_2', 'l3_0'
               , 'l4_0', 'l5_0', 'l5_1_2', 'l6_0', 'l6_1_1'
               , 'l7_0', 'l7_1_1', 'l8_0', 'l8_1_1', 'l9_0', 'l9_1_2', 'l10_0'
               , 'l10_1_2', 'l11_0', 'l11_1_2')
  expect_equal(kernel.parm.names$gammas, gammas.attr2)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.attr2.DCM)
  
  # DINA
  model.type <- "DINA"
  nattributes <- 3
  nclasses <- 2^nattributes
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  lambdas.DINA <- c('l1_0', 'l1_1_e', 'l2_0', 'l2_1_e', 'l3_0', 'l3_1_e', 'l4_0', 'l4_1_e'
                    , 'l5_0', 'l5_1_e', 'l6_0', 'l6_1_e', 'l7_0', 'l7_1_e', 'l8_0', 'l8_1_e'
                    , 'l9_0', 'l9_1_e', 'l10_0', 'l10_1_e', 'l11_0', 'l11_1_e')
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.DINA)
  # doesn't change with an added interaction
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.DINA)
  
  # DINO
  model.type <- "DINO"
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  lambdas.DINO <- c('l1_0', 'l1_1_e', 'l2_0', 'l2_1_e', 'l3_0', 'l3_1_e', 'l4_0', 'l4_1_e'
                    , 'l5_0', 'l5_1_e', 'l6_0', 'l6_1_e', 'l7_0', 'l7_1_e', 'l8_0', 'l8_1_e'
                    , 'l9_0', 'l9_1_e', 'l10_0', 'l10_1_e', 'l11_0', 'l11_1_e')
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.DINO)
  # doesn't change with an added interaction
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.DINO)
  
  # CRUM
  model.type <- "CRUM"
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  lambdas.CRUM <- c("l1_0", "l1_1_1", "l2_0", "l2_1_2", "l3_0", "l3_1_3", "l4_0", "l4_1_3"
                    , "l5_0", "l5_1_2", "l6_0", "l6_1_1", "l7_0", "l7_1_1", "l8_0", "l8_1_1"
                    , "l9_0", "l9_1_2", "l10_0", "l10_1_2", "l11_0", "l11_1_2")
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.CRUM)
  # Adding in an interaction, interaction term never appears, just all main effects
  lambdas.CRUM.interaction <- c("l1_0", "l1_1_1", "l2_0", "l2_1_2", "l2_1_1", "l3_0", "l3_1_3", "l4_0", "l4_1_3"
                    , "l5_0", "l5_1_2", "l6_0", "l6_1_1", "l7_0", "l7_1_1", "l8_0", "l8_1_1"
                    , "l9_0", "l9_1_2", "l10_0", "l10_1_2", "l11_0", "l11_1_2")
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.CRUM.interaction)
  
  # NIDO
  model.type <- "NIDO"
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  # main effects appear only once, as they are held constant across items
  lambdas.NIDO <- c("l1_0", "l1_1_1", "l2_0", "l1_1_2", "l3_0", "l1_1_3", "l4_0"
                    , "l5_0", "l6_0", "l7_0", "l8_0", "l9_0", "l10_0", "l11_0")
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.NIDO)
  # Adding in an interaction, interaction term never appears, just all main effects, thus unique lamdas are the same
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  expect_equal(kernel.parm.names$lambdas$lambda.names, lambdas.NIDO)
  
  # NCRUM
  model.type <- "NCRUM"
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test, nattributes, model.type)
  expect_is(kernel.parm.names, 'list')
  expect_equal(kernel.parm.names$gammas, gammas.attr3)
  pis.rias.NCRUM <- c("pi.1.star", "pi.2.star", "pi.3.star", "pi.4.star", "pi.5.star", "pi.6.star"
                      , "pi.7.star", "pi.8.star", "pi.9.star", "pi.10.star", "pi.11.star", "r.1.1.star"
                      , "r.2.2.star", "r.3.3.star", "r.4.3.star", "r.5.2.star", "r.6.1.star", "r.7.1.star"
                      , "r.8.1.star", "r.9.2.star", "r.10.2.star", "r.11.2.star")
  expect_equal(kernel.parm.names$pis.rias, pis.rias.NCRUM)
  # adding in interaction (now 2 rias for item with interaction)
  kernel.parm.names <- GetKernelParameterNames(qmatrix.test.interaction, nattributes, model.type)
  pis.rias.NCRUM.interaction <- c("pi.1.star", "pi.2.star", "pi.3.star", "pi.4.star", "pi.5.star"
                                  , "pi.6.star", "pi.7.star", "pi.8.star", "pi.9.star", "pi.10.star"
                                  , "pi.11.star", "r.1.1.star", "r.2.1.star", "r.2.2.star", "r.3.3.star"
                                  , "r.4.3.star", "r.5.2.star", "r.6.1.star", "r.7.1.star", "r.8.1.star"
                                  , "r.9.2.star", "r.10.2.star", "r.11.2.star")
  expect_equal(kernel.parm.names$pis.rias, pis.rias.NCRUM.interaction)
  
})