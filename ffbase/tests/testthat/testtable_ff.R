library(testthat)

context("table.ff")

test_that("table.ff works for one factor",{
   a <- factor(rep(c("A","B","C"), 10))
   af <- ff(a)

   t <- table(a)
   ft <- table.ff(af)
   
   names(dimnames(t)) <- NULL
   names(dimnames(ft)) <- NULL
   
   expect_equal(t, ft)
})  

test_that("table.ff works for two factors",{
   A <- factor(rep(c("A","B","C"), 10))
   a <- factor(rep(c("a","b"), 15))

   Af <- ff(A)
   af <- ff(a)

   t <- table(A,a)
   ft <- table.ff(Af,af)
   
   names(dimnames(t)) <- NULL
   names(dimnames(ft)) <- NULL
   
   expect_equal(t, ft)
})  

test_that("table.ff has same names as table.default",{
  iris_ffdf <- as.ffdf(iris)
  
  expected <- table(iris$Species)
  object <- table(iris_ffdf$Species)
  expect_equal(object, expected)
  expect_equal(dimnames(object), dimnames(expected))
  
  # only works with 'with' from 'with' branch
  #expect_equal(  with(iris_ffdf, table(Species))
  #             , with(iris, table(Species)))
})

# test_that("table.ff works for a large factor",{
#    x <- ff(length=1e7, vmode="integer")
#    for (i in chunk(x)){
#       x[i] <- repfromto(c(1,2,3,3,3), i[1], i[2])
#    }
#    levels(x) <- c("A","B","C") 
#    expect_equal(as.integer(table.ff(x)), c(2000000, 2000000, 6000000))
# })