context("r2_Wf")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

# loading bodyweight
data("bodyweight")

test_that("point estimate correct", {
    expect_equal(r2_Wf(msats, bodyweight, family = gaussian)$r2_Wf_full, 0.434, tolerance = 0.001)
    expect_equal(r2_Wf(snps, bodyweight, family = gaussian, type = "snps")$r2_Wf_full, 0.1285889, tolerance = 0.001)
})

test_that("matrix input works", {
    expect_equal(r2_Wf(as.matrix(msats), bodyweight, family = gaussian)$r2_Wf_full, 0.434, tolerance = 0.001)
    expect_equal(r2_Wf(as.matrix(snps), bodyweight, family = gaussian, type = "snps")$r2_Wf_full, 0.1285889, tolerance = 0.001)
})

test_that("bootstrapping works for msats", {
    expect_equal(r2_Wf(msats, trait = bodyweight, nboot = NULL,type = "msats")$r2_Wf_boot, NA)
    expect_equal(is.numeric(r2_Wf(msats, trait = bodyweight, nboot = 10,type = "msats")$r2_Wf_boot), TRUE)
    expect_equal(length(r2_Wf(msats, trait = bodyweight, nboot = 10,type = "msats")$r2_Wf_boot), 11)
})

test_that("bootstrapping works for snps", {
    expect_equal(r2_Wf(snps, trait = bodyweight, nboot = NULL,type = "snps")$r2_Wf_boot, NA)
    expect_equal(is.numeric(r2_Wf(snps, trait = bodyweight, nboot = 10,type = "snps")$r2_Wf_boot), TRUE)
    expect_equal(length(r2_Wf(snps, trait = bodyweight, nboot = 10,type = "snps")$r2_Wf_boot), 11)
})