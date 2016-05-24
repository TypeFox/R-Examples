context("r2_hf")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

test_that("nboot = 0 works", {
    expect_equal(r2_hf(msats, type = "msats")$r2_hf_full, 0.2797, 
                 tolerance = 0.001)
    expect_equal(r2_hf(snps, type = "snps")$r2_hf_full, 0.6698, 
                 tolerance = 0.001)
})

# test_that("subsets = NULL works", {
#     expect_equal(r2_hf(msats, nboot_loci = 0, type = "msats")$r2_hf_full, 0.2797, tolerance = 0.001)
#     expect_equal(r2_hf(snps, nboot_loci = 0, type = "snps")$r2_hf_full, 0.6698, tolerance = 0.001)
#     expect_equal(r2_hf(msats, nboot_loci = 5, type = "msats")$r2_hf_full, 0.2797, tolerance = 0.001)
#     expect_equal(r2_hf(snps, nboot_loci = 5, type = "snps")$r2_hf_full, 0.6698, tolerance = 0.001)
# })

test_that("Matrix input works", {
    expect_equal(r2_hf(as.matrix(msats),type = "msats")$r2_hf_full, 
                 0.2797, tolerance = 0.001)
    expect_equal(r2_hf(as.matrix(snps), type = "snps")$r2_hf_full, 
                 0.6698, tolerance = 0.001)
})

test_that("bootstrapping works for msats", {
    expect_equal(r2_hf(msats, nboot = NULL,type = "msats")$r2_hf_boot, NA)
    expect_equal(is.numeric(r2_hf(msats, nboot = 10,type = "msats")$r2_hf_boot), TRUE)
    expect_equal(length(r2_hf(msats, nboot = 10,type = "msats")$r2_hf_boot), 11)
})

test_that("bootstrapping works for snps", {
    expect_equal(r2_hf(snps, nboot = NULL,type = "snps")$r2_hf_boot, NA)
    expect_equal(is.numeric(r2_hf(snps, nboot = 10,type = "snps")$r2_hf_boot), TRUE)
    expect_equal(length(r2_hf(snps, nboot = 10,type = "snps")$r2_hf_boot), 11)
})

