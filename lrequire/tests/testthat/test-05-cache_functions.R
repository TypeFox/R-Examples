context('cache_functions')

tmp.file <- 'dyntest_cache_functions.R'
cache <- get.module.cache()

create.tmp.R.module <- function(exports.value, change_code) {
  fileConn <- file(tmp.file)
  writeLines(c('value <- 1'), fileConn)
  close(fileConn)
}

test_that('reset.module.cache removes appropriate cache entries', {
  reset.module.cache()

  # At this point there should be two hidden items in the cache:
  #  .:module.cache
  #  .:module.change_code
  #  .:warn.not.found
  expect_equal(length(cache), 3)
  expect_equal(length(ls(cache)), 0)
})

test_that('when lrequire-ing a file, there should be two additional items in the cache', {
  reset.module.cache()
  create.tmp.R.module()

  lrequire(tmp.file, character.only = TRUE)
  expect_equal(length(ls(cache)), 2)
})

test_that('remove.from.module.cache should remove the file from the cache', {
  remove.from.module.cache(tmp.file, character.only = TRUE)

  expect_equal(length(ls(cache)), 0)
})

# cleanup temporary file
file.remove(tmp.file)
