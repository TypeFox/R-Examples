context("MLH")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

test_that("sMLH is correct", {
    expect_equal(mean(MLH(msats)), 0.5959596, tolerance = 0.001)
    expect_equal(sd(MLH(msats)), 0.1663811, tolerance = 0.001)
})