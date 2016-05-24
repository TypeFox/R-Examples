require(testthat)


context("Parsing")

test_that("parseEdits all works",{
  x <- c( "2*x < 1"
        , "if (A=='a') B == 'b'"
        , "if (A=='a') B == FALSE"
        , "if (A=='a') B > 1"
        , "if (c==1) B || C == FALSE"
        )
  e <- parseEdits(x)  
  expect_equal(length(e), 5)
})

test_that("parseEdits num works",{
  x <- c( "2*x < 1"
        , "if (A=='a') B == 'b'"
        , "if (A=='a') B == FALSE"
        , "if (A=='a') B > 1"
        , "if (c==1) B || C == FALSE"
        )
  e <- parseEdits(x, "num")  
  expect_equal(length(e), 1)
  u <- new.env(); u$x = 1
  expect_equal(eval(e,u),eval(expression(2*x < 1), u))
})

test_that("parseEdits cat works",{
  x <- c( "2*x < 1"
        , "if (A=='a') B == 'b'"
        , "if (A=='a') B == FALSE"
        , "if (A=='a') B > 1"
        , "if (c==1) B || C == FALSE"
        )
  e <- parseEdits(x, "cat")  
  expect_equal(length(e), 2)
# test fails: equivalence of expressions cannot be tested like this. 
#  expect_equivalent(e, expression( if (A == "a") B == "b"
#                            , if (A == "a") B == FALSE
#                            )
#              )
})

test_that("parseEdits mix works",{
  x <- c( "2*x < 1"
        , "if (A=='a') B == 'b'"
        , "if (A=='a') B == FALSE"
        , "if (A=='a') B > 1"
        , "if (c==1) B || C == FALSE"
        )
  e <- parseEdits(x, "mix")  
  expect_equal(length(e), 2)
# test fails: equivalence of expressions cannot be tested like this. 
#  expect_equivalent(e, expression( if (A == "a") B > 1
#                            , if (c == 1) B || C == FALSE
#                            )
#               , label=deparse(e)
#               )
})

context("editfile")
test_that("assignments are parsed by editfile",{
  e <- editfile("edit_test_1.txt",type='num')
  expect_equivalent(getA(e),array(c(-10,1),dim=c(1,2)))
})

test_that("empty files are parsed correctly",{
  expect_equal(
    nedits(editfile(textConnection("#test"),type='num'))    
    , 0
  )
})



