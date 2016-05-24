context('reload_and_caching')

tmp.file <- 'dyntest_reload_and_caching.R'

create.tmp.R.module <- function(exports.value, change_code) {
  fileConn <- file(tmp.file)
  writeLines(c(paste("module.exports <- ", exports.value), paste("module.change_code = ", change_code)), fileConn)
  close(fileConn)
}

test_that('when seen for the first time, it should load file from cache', {
  create.tmp.R.module(exports.value = 1, change_code = 0)
  ret.val <- lrequire(dyntest_reload_and_caching)

  expect_equal(ret.val, 1)
})

test_that('when seen previously, it should load file from cache', {
  create.tmp.R.module(exports.value = 2, change_code = 1)
  ret.val <- lrequire('dyntest_reload_and_caching.R')

  expect_equal(ret.val, 1)
})

test_that('when forcing a reload, it should load from disk', {
  create.tmp.R.module(exports.value = 3, change_code = 1)
  ret.val <- lrequire(file.path(getwd(), './././dyntest_reload_and_caching'), force.reload = TRUE)

  expect_equal(ret.val, 3)
})

test_that('since change_code = 1 and mtime has changed, it should load from disk', {
  Sys.sleep(1)  # Need to delay so the modification time is different than the previous file
  create.tmp.R.module(exports.value = 4, change_code = 1)
  ret.val <- lrequire('./././dyntest_reload_and_caching')

  expect_equal(ret.val, 4)
})

# remove temporary file
file.remove(tmp.file)
