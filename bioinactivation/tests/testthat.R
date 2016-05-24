
library(testthat)
library(bioinactivation)

#==============================================================================

#- TEST PREDICTIONS FOR ISOTHERMAL PROFILES (WITH KNOWN SOLUTIONS)

test_that("Prediction Bigelow", {

    temp_profile <- data.frame(time = c(0, 5), temperature = c(100, 100))

    simulation_model <- "Bigelow"
    times <- seq(0, 2, length = 100)
    parms <- list(D_R = 5, z = 10, temp_ref = 90, N0 = 1e5)

    bigelow_results <- predict_inactivation(simulation_model, times,
                                            parms, temp_profile)


    expect_equal(bigelow_results$simulation$logN[1], 5, tolerance = 1e-4)
    expect_equal(bigelow_results$simulation$logS[100], -4, tolerance = 1e-4)
})

#------------------------------------------------------------------------------

test_that("Prediction Peleg", {

    temp_profile <- data.frame(time = c(0, 5), temperature = c(100, 100))

    simulation_model <- "Peleg"
    times <- seq(0, 2, length = 100)
    parms <- list(k_b = 0.1, temp_crit = 95, n = 2, N0 = 1e5)

    peleg_results <- predict_inactivation(simulation_model, times,
                                          parms, temp_profile)

    expect_equal(peleg_results$simulation$logN[1], 5, tolerance = 1e-4)
    expect_equal(peleg_results$simulation$logS[100], -3.9, tolerance = 1e-4)
})

#------------------------------------------------------------------------------

test_that("Prediction Mafart", {

    temp_profile <- data.frame(time = c(0, 5), temperature = c(100, 100))

    simulation_model <- "Mafart"
    times <- seq(0, 2, length = 100)
    parms <- list(delta_ref = 10, temp_ref = 90, z = 10, p = 2, N0 = 1e5)

    mafart_results <- predict_inactivation(simulation_model, times,
                                           parms, temp_profile)

    expect_equal(mafart_results$simulation$logN[1], 5, tolerance = 1e-4)
    expect_equal(mafart_results$simulation$logS[100], -4, tolerance = 1e-4)
})

#==============================================================================

