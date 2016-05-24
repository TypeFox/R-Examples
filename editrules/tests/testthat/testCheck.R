library(testthat)

context("Edit checking")

test_that("Checking a data set works",{
   edt <- editmatrix("2*x==y")
   dat <- data.frame(x=c(2,1), y=c(4,5))
   cr <- checkRows(edt, dat)
   expect_equal(cr, c(TRUE,FALSE))
})

test_that("Checking a vector works",{
   edt <- editmatrix(editrules = "2*x==y")
   dat <- c(x=2, y=1)
   expect_equal(as.logical(violatedEdits(edt, dat)), c(TRUE))
})

test_that("Showing data errors works",{
   edt <- editmatrix(editrules = "2*x==y\n5*x==y")
   dat <- data.frame(x=c(2,1), y=c(4,5))
   errors <- violatedEdits(edt, dat)
   
   # rule names equal?
   expect_equal(colnames(errors), rownames(edt))
   
   dimnames(errors) <- NULL
   expect_equal(errors[,,drop=FALSE], matrix( c( FALSE, TRUE
                                 , TRUE , FALSE
                                 )
                              , nrow=2
                              , byrow=TRUE
                              
                              )
               )   
})


test_that("Error lists works",{
   edt <- editmatrix(editrules = "2*x==y\n5*x==y")
   dat <- data.frame(x=c(2,1), y=c(4,5))
   errors <- listViolatedEdits(edt, dat)
   #expect_equal(errors, list("1"=c(2), "2"=c(1)))
   #print(str(errors))
   #print(showErrors(dat, edt)) 
})
