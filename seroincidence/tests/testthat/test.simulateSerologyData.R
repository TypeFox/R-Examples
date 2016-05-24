context("Simulation")

test_that("simulateSerologyData produces samples", {
    expect_that(simulateSerologyData(n = 100), is_a("data.frame"))
    expect_that(nrow(simulateSerologyData(n = 100)), equals(100))
})


