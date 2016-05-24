context("chunk")

test_that("chunk", {
  # normal chunk.size
  x = 1:9
  ch = chunk(x, chunk.size=3)
  expect_equal(ch, list(1:3, 4:6, 7:9))

  # normal n.chunks
  x = 1:9
  ch = chunk(x, n.chunks=3)
  expect_equal(ch, list(1:3, 4:6, 7:9))

  # chunk.size uneven
  x = 1:10
  ch = chunk(x, chunk.size=3)
  expect_equal(ch, list(1:3, 4:6, 7:8, 9:10))

  # n.chunks uneven
  ch = chunk(1:9, n.chunks=4)
  expect_equal(length(ch), 4)

  x = letters[1:10]
  ch = chunk(x, n.chunks = 2)
  expect_equal(ch, list(letters[1:5], letters[6:10]))

  # errors
  x = letters[1:10]
  expect_error(chunk(x, chunk.size=1, n.chunks=3))
  expect_error(chunk(x, chunk.size=1:2))
  expect_error(chunk(x, n.chunks=list()))

  x = as.list(letters[1:10])
  ch = chunk(x, chunk.size=5)
  expect_equal(ch, list(as.list(letters[1:5]), as.list(letters[6:10])))

  x = letters
  ch = chunk(x, chunk.size=4, shuffle=TRUE)
  expect_equal(sort(letters), sort(unlist(ch)))
  expect_true(all(sapply(ch, length) %in% c(3, 4)))

  # test that smaller levels get chosen randomly
  x = 1:5
  counts = sapply(1:100, function(i) {
    ch = chunk(x, chunk.size=3, shuffle=TRUE)
    sapply(ch, length) == 2
  })
  counts = rowSums(counts)
  expect_true(all(counts > 30))

  # test proportions
  x = 1:10
  ch = chunk(x, props = c(3, 7))
  expect_equal(sapply(ch, length), c(3, 7))
  expect_equal(unlist(ch), x)

  expect_true(length(chunk(x, props=1)) == 1L)

})

