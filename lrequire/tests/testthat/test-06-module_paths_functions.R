context('cache_functions')

orig.paths <- get.module.paths()

test_that('without altering module.paths, get.module.paths returns the same vector', {
  paths <- get.module.paths()

  expect_equal(paths, orig.paths)
})

test_that('appending path to end of module.paths should function as expected', {
  append.module.paths("/")

  paths <- get.module.paths()
  expect_equal(length(paths), length(orig.paths)+1)

  expect_equal(paths[length(orig.paths)+1], "/")
})

test_that('removal of previously added path should function as expected', {
  remove.module.paths(length(orig.paths)+1)

  paths <- get.module.paths()
  expect_equal(length(paths), length(orig.paths))

  expect_equal(paths, orig.paths)
})

test_that('inserting path into front of module.paths should function as expected', {
  append.module.paths("/", 0)

  paths <- get.module.paths()
  expect_equal(length(paths), length(orig.paths)+1)

  expect_equal(paths[1], "/")
})

test_that('inserting path into middle of module.paths should function as expected', {
  append.module.paths("/", 3)

  paths <- get.module.paths()
  expect_equal(length(paths), length(orig.paths)+2)

  expect_equal(paths[4], "/")
})

test_that('removal of previously added paths should function as expected', {
  remove.module.paths(1, 4)

  paths <- get.module.paths()
  expect_equal(length(paths), length(orig.paths))

  expect_equal(paths, orig.paths)
})
