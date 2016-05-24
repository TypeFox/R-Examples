
test_that( "not in operator works", {
  expect_true( ! "aa" %!in% c("aa", "bb"))
  expect_true( "aa" %!in% c("bb","cc"))
})

test_that("%without% works", {
  expect_equal( c("A","B") %without% "A", "B" )
  expect_equal( c("A","B") %without% "C", c("A","B") )
})

test_that("%of% works", {
  expect_true( iris %of% "data.frame")
  expect_true( glm(Sepal.Length~Species, data = iris) %of% "glm")
  expect_true( glm(Sepal.Length~Species, data = iris) %of% "lm")
  expect_true( glm %of% "function")
  expect_true( (x~y) %of% "formula")
})


