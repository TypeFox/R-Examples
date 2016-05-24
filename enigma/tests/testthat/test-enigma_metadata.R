# tests for enigma_metadata fxn in taxize
context("enigma_metadata")

## Domestic Market Flight Statistics (Final Destination)
# enigma_metadata('us.gov.dot.rita.trans-stats.air-carrier-statistics.t100d-market-all-carrier')

test_that("enigma_metadata basic functionality works", {
  skip_on_cran()
  
  res1 <- enigma_metadata(dataset='us.gov.dot.rita.trans-stats.air-carrier-statistics')
  expect_is(res1, "enigma_meta")
  expect_true(res1$success)
  expect_is(res1$datapath, "character")
  expect_is(res1$info, "list")
  expect_is(res1$info$paths, "list")
})

test_that("enigma_metadata parent node data differs from child node data", {
  skip_on_cran()
  
  parent <- enigma_metadata(dataset='us.gov.whitehouse')
  child <- enigma_metadata(dataset='us.gov.whitehouse.visitor-list')
  expect_equal(names(child$info), c('info','table','ancestor_datapaths','db_boundary_datapath','db_boundary_label','db_boundary_tables','paths'))
  expect_equal(names(parent$info), c("paths","immediate_nodes","children_tables"))
})
