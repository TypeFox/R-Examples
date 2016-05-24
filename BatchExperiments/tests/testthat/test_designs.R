context("designs")

test_that("designs", {
  # one point, empty function call
  d = BatchExperiments:::designIterator(ex=list())
  expect_true(d$hasNext())
  expect_equal(d$n.states, 1)
  expect_equal(d$nextElem(), list())

  d = BatchExperiments:::designIterator(ex=list(), .design=data.frame(a=1:3))
  expect_true(d$hasNext())
  expect_equal(d$n.states, 3L)
  xs = replicate(3, d$nextElem(), simplify=FALSE)
  expect_equal(xs, list(list(a=1L), list(a=2L), list(a=3L)))

  d = BatchExperiments:::designIterator(ex=list(), .design=data.frame(a=1:2, b=4:5))
  expect_true(d$hasNext())
  expect_equal(d$n.states, 2L)
  xs = replicate(2, d$nextElem(), simplify=FALSE)
  expect_equal(xs, list(list(a=1L, b=4L), list(a=2L, b=5L)))

  d = BatchExperiments:::designIterator(ex=list(a=c("x", "y"), b=4:5), .design=data.frame(c=c(1,6), d=c(7,9)))
  expect_true(d$hasNext())
  expect_equal(d$n.states, 8L)
  xs = replicate(8, as.data.frame(d$nextElem()), simplify=FALSE)
  xs = do.call(rbind, xs)
  data = data.frame(
    c = c(1, 1, 1, 1, 6, 6, 6, 6),
    d = c(7, 7, 7, 7, 9, 9, 9, 9),
    a = c("x", "y", "x", "y", "x", "y", "x", "y"),
    b = c(4, 4, 5, 5, 4, 4, 5, 5)
  )
  expect_equal(xs, data[names(xs)])
  expect_false(d$hasNext())
  d$reset()
  expect_equal(d$n.states, 8L)
  xs = replicate(8, as.data.frame(d$nextElem()), simplify=FALSE)
  xs = do.call(rbind, xs)
  expect_equal(xs, data[names(xs)])

  d = BatchExperiments:::designIterator(ex=list(a=as.factor(c("u", "v"))), .design=data.frame(b=as.factor(c("a", "b")), c=c("x", "y"), stringsAsFactors=FALSE))
  x  = d$nextElem()
  y = list(b=factor("a", levels=c("a", "b")), c="x", a=factor("u", levels=c("u", "v")))
  expect_equal(x, y[names(x)])

  #expect_equal(nrow(d$design), 2)
  #expect_equal(ncol(d$design), 1)
  #expect_equal(colnames(d$design), c("a"))
  #expect_true(is.numeric(d$design$a))
  #d = makeDesign(a, design=data.frame(a=1:2, b=c("a", "b")))
  #expect_equal(nrow(d$design), 2)
  #expect_equal(ncol(d$design), 2)
  #expect_equal(colnames(d$design), c("a", "b"))
  #expect_true(is.numeric(d$design$a) && is.character(d$design$b))
  #d = makeDesign(a, design=data.frame(a=1:2, b=c("a", "b")), exhaustive=list(c=pi, d=1:5))
  #expect_equal(nrow(d$design), 10)
  #expect_equal(ncol(d$design), 4)
  #expect_equal(colnames(d$design), c("a", "b", "c", "d"))
  #expect_true(is.numeric(d$design$a) && is.character(d$design$b)
  #  && is.numeric(d$design$c) && is.numeric(d$design$d))

  #d = check_error(makeDesign(a, exhaustive=list(x=iris)), "primitive")
})
