context("make_taxon_fromclass")

test_that("make_taxon_fromclass works", {
  df <- data.frame(rank=c('family','tribe','subtribe','genus','subgenus','species'),
      name=c('Helianthi','Helianthi','Helianthi','Poa','Festuci','Poa annua'),
      id=c(1,2,3,4,5,6), stringsAsFactors = FALSE)
  aa <- apply(df, 1, make_taxon_fromclass)

  expect_is(aa, "list")
  expect_is(aa[[1]], "list")
  expect_is(aa[[1]]$family, "taxonref")
  expect_named(aa[[2]]$tribe, c('rank', 'name', 'id', 'uri'))
})
