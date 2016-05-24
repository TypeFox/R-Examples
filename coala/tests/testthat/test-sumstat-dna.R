context("SumStat DNA")

test_that("DNA can be simulated", {
  if (!has_seqgen()) skip("seq-gen not installed")
  model <- coal_model(c(5, 5), 1, 10) +
    locus_single(5) +
    feat_pop_merge(.5, 2, 1) +
    feat_mutation(5, model = "GTR", gtr_rates = 1:6) +
    feat_recombination(par_const(1)) +
    sumstat_dna("dna")

  stat <- simulate(model)
  expect_that(stat$dna, is_a("list"))
  expect_equal(length(stat$dna), 2)
  expect_is(stat$dna[[1]], "matrix")
  expect_is(stat$dna[[2]], "matrix")

  temp_files_before <- list.files(tempdir(), pattern = "^coala-[0-9]+-")
  expect_error(simulate(model + sumstat_sfs() + feat_outgroup(2)))
  expect_error(simulate(model + locus_trio()))

  # Remove tempfiles that may remain because of error exit
  temp_files_after <- list.files(tempdir(), pattern = "^coala-[0-9]+-")
  temp_files_diff <- temp_files_after[!temp_files_after %in% temp_files_before]
  unlink(file.path(tempdir(), temp_files_diff))
})
