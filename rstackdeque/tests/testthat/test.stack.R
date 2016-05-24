test_that("peek_top assignment works as expected", {
  s <- rstack()
  s <- insert_top(s, data.frame(a = 1, b = 2))
  s <- insert_top(s, data.frame(a = 3, b = 4))
  peek_top(s)$a <- 300
  
  slist <- as.list(s)
  expect_that(slist[[1]], equals(data.frame(a = 300, b = 4)))
  
  peek_top(s) <- data.frame(a = 300, b = 400)
  slist <- as.list(s)
  expect_that(slist[[1]], equals(data.frame(a = 300, b = 400)))
})


test_that("insert and remove work as expected", {
  s <- rstack()
  s <- insert_top(s, "c")
  s <- insert_top(s, "b")
  s <- insert_top(s, "a")
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("a", "b", "c"))) )
  
  ## drop one
  sp <- without_top(s)
  expect_that(length(sp), equals(2))
  expect_that(as.list(sp), is_identical_to(as.list(c("b", "c"))) )  
  
  ## s still ok?
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("a", "b", "c"))) )
  
  ## add one
  se <- insert_top(s, "aa")
  expect_that(length(se), equals(4))
  expect_that(as.list(se), is_identical_to(as.list(c("aa", "a", "b", "c"))) )  
  
  ## s still ok?
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("a", "b", "c"))) )
})