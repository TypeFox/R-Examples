test_that("peek_front assignment works as expected", {
  d <- rdeque()
  d <- insert_front(d, data.frame(a = 1, b = 2))
  d <- insert_front(d, data.frame(a = 3, b = 4))
  peek_front(d)$a <- 300
  
  dlist <- as.list(d)
  expect_that(dlist[[1]], equals(data.frame(a = 300, b = 4)))
  
  peek_front(d) <- data.frame(a = 300, b = 400)
  dlist <- as.list(d)
  expect_that(dlist[[1]], equals(data.frame(a = 300, b = 400)))
})

test_that("peek_back assignment works as expected", {
  d <- rdeque()
  d <- insert_front(d, data.frame(a = 1, b = 2))
  d <- insert_front(d, data.frame(a = 3, b = 4))
  peek_back(d)$a <- 100
  
  dlist <- as.list(d)
  expect_that(dlist[[2]], equals(data.frame(a = 100, b = 2)))
  
  peek_back(d) <- data.frame(a = 100, b = 2)
  dlist <- as.list(d)
  expect_that(dlist[[2]], equals(data.frame(a = 100, b = 2)))
})


test_that("insert_front and without_front work as expected", {
  s <- rdeque()
  s <- insert_front(s, "c")
  s <- insert_front(s, "b")
  s <- insert_front(s, "a")
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("a", "b", "c"))) )
  
  ## drop one
  sp <- without_front(s)
  expect_that(length(sp), equals(2))
  expect_that(as.list(sp), is_identical_to(as.list(c("b", "c"))) )  
  
  ## s still ok?
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("a", "b", "c"))) )
  
  ## add one
  se <- insert_front(s, "aa")
  expect_that(length(se), equals(4))
  expect_that(as.list(se), is_identical_to(as.list(c("aa", "a", "b", "c"))) )  
  
  ## s still ok?
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("a", "b", "c"))) )

  s <- without_front(s)
  s <- without_front(s)
  s <- without_front(s)
  expect_that(empty(s), equals(TRUE))
})


test_that("peek_front and peek_back work for nearly-empty deques", {
  d <- rdeque()
  d <- insert_front(d, "a")
  expect_that(peek_front(d), equals("a"))
  expect_that(peek_back(d), equals("a"))

  d2 <- rdeque()
  d2 <- insert_front(d2, "a")
  expect_that(peek_front(d2), equals("a"))
  expect_that(peek_back(d2), equals("a"))
  
  d3 <- rdeque()
  d3 <- insert_front(d3, "a")
  d3 <- insert_front(d3, "b")
  expect_that(peek_front(d3), equals("b"))
  expect_that(peek_back(d3), equals("a"))

  d4 <- rdeque()
  d4 <- insert_back(d4, "a")
  d4 <- insert_back(d4, "b")
  expect_that(peek_back(d4), equals("b"))
  expect_that(peek_front(d4), equals("a"))
})

test_that("insert_back and without_back work as expected", {
  s <- rdeque()
  s <- insert_back(s, "c")
  s <- insert_back(s, "b")
  s <- insert_back(s, "a")
  expect_that(length(s), equals(3))
  expect_that(as.list(s), is_identical_to(as.list(c("c", "b", "a"))) )
  
  ## drop one
  sp <- without_back(s)
  expect_that(length(sp), equals(2))
  expect_that(as.list(sp), is_identical_to(as.list(c("c", "b"))) )  
  
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
  
  s <- without_back(s)
  s <- without_back(s)
  s <- without_back(s)
  expect_that(empty(s), equals(TRUE))
  
})
