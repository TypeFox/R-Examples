context("Calculating identity score")

# test_that("rpp_align() produces correct results", {
#   prepared_ref <- readRDS("prepared_ref.rda")
#   prepared_com <- readRDS("prepared_com.rda")
#   results <- readRDS("results.rda")
#   means <- readRDS("means.rda")
# 
#   scores <- rcpp_align(prepared_ref,prepared_com)
# 
#   expect_equal(class(scores),"list")
#   sresults = scores$results
#   expect_equal(sresults[1,],results[1,])
#   smeans = scores$means
#   expect_equal(smeans,means)
# 
# })


test_that("compare_alignments() produces a list", {
  data(reference_alignment)
  data(comparison_alignment)

  scores <- compare_alignments(reference_alignment,comparison_alignment)

  expect_equal(class(scores),"list")
})


# test_that("prepare_alignment_matrix() produces correct outputs",{
#   data(reference_alignment)
#   prepared_ref <- readRDS("prepared_ref.rda")
#   r2 <- prepare_alignment_matrix(reference_alignment)
# 
#   expect_equal(prepared_ref,r2)
# 
# })



