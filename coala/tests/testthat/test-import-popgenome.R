context("Import PopGenome")

test_that("importing PopGenome Data works", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()

  seg_sites <- as.segsites(data_pg)
  expect_equal(length(seg_sites), 1)
  expect_true(is_segsites(seg_sites[[1]]))
  expect_equal(nrow(seg_sites[[1]]), 12)

  expect_false(is.null(get_positions(seg_sites[[1]])))
  expect_true(all(get_positions(seg_sites[[1]]) >= 0))
  expect_true(all(get_positions(seg_sites[[1]]) <= 1))
})
