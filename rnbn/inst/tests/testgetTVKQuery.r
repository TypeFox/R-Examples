context("Test getTVKQuery")
dt <- getTVKQuery(query='Badger')
test_that("query works1", {
    expect_that('NHMSYS0000080191' %in% dt$ptaxonVersionKey, is_true()) 
})
dt <- getTVKQuery(query='myotis')
test_that("query works2", {
    expect_that(!'Genus' %in% dt$rank, is_true())
    expect_that('Synonym' %in% dt$nameStatus, is_true())
})
dt <- getTVKQuery(query='myotis',species_only=F,rec_only=T)
test_that("query works2", {
    expect_that('Genus' %in% dt$rank, is_true())
    expect_that(!'Synonym' %in% dt$nameStatus, is_true())
})
dt <- getTVKQuery(query='myotis', top = T)
test_that("query works2", {
    expect_that(nrow(dt) == 1, is_true())
    expect_that(dt$searchMatchTitle == 'Myotis myotis', is_true())
})