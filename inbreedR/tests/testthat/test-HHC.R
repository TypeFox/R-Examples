context("HHC")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

test_that("HHC mean is in range 0 - 1", {
    expect_equal(mean(round(HHC(msats, reps = 100, CI = 0.95)$HHC_vals)), 0.5, tolerance = 0.5)
    # expect_equal(mean(round(HHC(snps, niter = 100, CI = 0.95)$HHC_vals)), 0.5, tolerance = 0.5)
})

