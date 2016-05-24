# skip_on_os(windows)
library(futureheatwaves)
context("Identifying heatwaves")

data("datafr", package = "futureheatwaves")
hws_r <- IDHeatwavesR(datafr = datafr, threshold = 84.632, numDays = 2)
hws_cpp <- IDHeatwavesCPPwrapper(datafr = datafr, threshold = 84.632,
                                 numDays = 2)

test_that("C++ and R ID heatwaves identify same heatwaves", {
        expect_equal(hws_r$hw, hws_cpp$hw)
        expect_equal(hws_r$hw.number, hws_cpp$hw.number)
})

hw_data <- data.frame(date = rep("A", 8),
                      tmpd = c(90, 95, 98, 94.9, 100, 95, 101, 92))
hw <- c(0, 1, 1, 0, 1, 1, 1, 0)
hw.number <- c(0, 1, 1, 0, 2, 2, 2, 0)
R_fun_ids <- IDHeatwavesR(datafr = hw_data, threshold = 95, numDays = 2)
Cpp_fun_ids <- IDHeatwavesCPPwrapper(datafr = hw_data, threshold = 95,
                                     numDays = 2)

test_that("Heat wave ID functions (R and C++) giving correct answers", {
        expect_equal(R_fun_ids$hw, hw)
        expect_equal(R_fun_ids$hw.number, hw.number)
        expect_equal(Cpp_fun_ids$hw, hw)
        expect_equal(Cpp_fun_ids$hw.number, hw.number)
})
