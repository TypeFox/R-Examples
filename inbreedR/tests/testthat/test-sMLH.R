context("sMLH")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

test_that("sMLH is correct", {
    expect_equal(mean(sMLH(msats)), 1.000855, tolerance = 0.001)
    expect_equal(sd(sMLH(msats)), 0.2791373, tolerance = 0.001)
})
