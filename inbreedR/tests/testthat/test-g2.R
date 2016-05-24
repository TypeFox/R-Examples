context("g2_functions")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]



test_that("g2 point estimates are correct", {
    expect_equal(round(g2_microsats(msats, nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02179805)
    expect_equal(round(g2_snps(snps, nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02270884)
    })

test_that("bootstrapping works", {
    expect_equal(sum(is.na(g2_microsats(msats, nperm = 0, nboot = 20, CI = 0.95)$CI_boot)), 0) # CI calculated ?
    expect_equal(sum(is.na(g2_snps(snps, nperm = 0, nboot = 20, CI = 0.95)$CI_boot)), 0)
})

test_that("permutation works", {
    expect_equal(sum(is.na(g2_microsats(msats, nperm = 20, nboot = 0, CI = 0.95)$p_val)), 0) # p_val calculated ?
    expect_equal(sum(is.na(g2_snps(snps, nperm = 20, nboot = 0, CI = 0.95)$p_val)), 0)
})

# test_that("parallelization works", {
#     expect_that(g2_snps(snps, nperm = 20, nboot = 0, CI = 0.95, parallel = TRUE,
#                         ncores = 2), 
#                 not(throws_error()))
# })

test_that("matrix input works", {
    expect_equal(round(g2_microsats(as.matrix(msats), nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02179805)
    expect_equal(round(g2_snps(as.matrix(snps), nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02270884)
})

# different missing values

msats_miss_num <- msats
msats_miss_num[is.na(msats_miss_num)] <- -99

snps_miss_num <- snps
snps_miss_num[is.na(snps_miss_num)] <- -99


test_that("g2 point estimates are correct", {
    expect_equal(round(g2_microsats(msats_miss_num, nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02179805)
    expect_equal(round(g2_snps(snps_miss_num, nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02270884)
})





