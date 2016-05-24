context("treebase")

test_that("search treebase works as expected", {
 skip_on_cran()         
 trees <- search_treebase("Huelsenbeck", "author")
 expect_that(is(trees, "list"), is_true())
 expect_that(is(trees[[1]], "phylo"), is_true())
})
