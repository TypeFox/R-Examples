
context("Deterministic corrections")
library(editrules)
.onLoad(0,0)
test_that("parser",{
   e <- correctionRules(expression(
      if ( x > y ) x <- y,
      if ( is.na(z) ) x <- 0
   ))
   expect_true(all(sapply(e,is.language)))
   expect_equal(length(e),2)
   expect_equal(sort(getVars(e)),c('x','y','z'))
   expect_error(correctionRules("if (is.na(x)) x <- mean(y)"))
   expect_equal(
      as.character(correctionRules(expression(if ( x < 0 ) x <- 0))),
      as.character(correctionRules("if ( x < 0 ) x <- 0",file=FALSE))
   )
   e <- correctionRules("if ( x < 0 ) x <- 0",file=FALSE)
   expect_equal(as.character(e), as.character(correctionRules(as.character(e),file=FALSE)))

})



test_that("correction",{
   expect_equal(
      correctWithRules(
         correctionRules("if ( x < 0 ) x <- 0",file=FALSE),
         data.frame(x=c(-1,0,1))
      )$corrected,
      data.frame(x=c(0,0,1))
   )
       
})

test_that("logging",{
   b <- correctionRules("if ( x < 0 ) x <- 0", file=FALSE)
   expect_equivalent(
      correctWithRules(
         b,
         data.frame(x=-1,0,1)
      )$corrections
    , data.frame(row=1,variable='x',old=format(-1),new=format(0),how=as.character(b),stringsAsFactors=FALSE)
   )
})


test_that("logging for changing a variable twice",{
   df <- data.frame(
      x = 1,
      y = NA,
      z = 2)

   u <- correctionRules(expression(
      if ( x == 1 ) z <- NA,
      if ( is.na(z) ) z <- 1
   ))
   x <- correctWithRules(u,df)

   expect_equal(x$corrections[2,'old'],format(NA))
   expect_equal(x$corrections[2,'new'],format(1))
})


test_that("changing multiple variables in a record with different rules",{

   df <- data.frame(
    x = 1,
   y = NA,
   z = 2)


   w <- correctionRules(expression(
      if ( x == 1 ) x <- NA,
      if ( is.na(y) ) y <- 1
   ))
   x <- correctWithRules(w,df)
   expect_equal(nrow(x$corrections),2)
   expect_equal(x$corrections$old,sapply(c(1,NA),format))
   expect_equal(x$corrections$new,sapply(c(NA,1),format))

})

