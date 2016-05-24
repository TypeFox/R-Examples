test_that("peek_front assignment works as expected", {
  d <- rpqueue()
  d <- insert_back(d, data.frame(a = 1, b = 2))
  d <- insert_back(d, data.frame(a = 3, b = 4))
  peek_front(d)$a <- 100
  
  dlist <- as.list(d)
  expect_that(dlist[[1]], equals(data.frame(a = 100, b = 2)))
  
  peek_front(d) <- data.frame(a = 100, b = 200)
  dlist <- as.list(d)
  expect_that(dlist[[1]], equals(data.frame(a = 100, b = 200)))
})


test_that("insert_back and without_front work as expected", {
  s <- rpqueue()
  s <- insert_back(s, "c")
  s <- insert_back(s, "b")
  s <- insert_back(s, "a")
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("c", "b", "a"))) )
  
  ## drop one
  sp <- without_front(s)
  expect_that(length(sp), equals(2))
  expect_that(as.list(sp), is_identical_to(as.list(c("b", "a"))) )  
  
  ## s still ok?
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("c", "b", "a"))) )
  
  ## add one
  se <- insert_back(s, "aa")
  expect_that(length(se), equals(4))
  expect_that(as.list(se), is_identical_to(as.list(c("c", "b", "a", "aa"))) )  
  
  ## s still ok?
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("c", "b", "a"))) )
  
  s <- without_front(s)
  s <- without_front(s)
  s <- without_front(s)
  expect_that(empty(s), equals(TRUE))
  
})

