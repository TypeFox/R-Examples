library(testthat)

context("Editmatrix attributes")

test_that("editrules can derive the correct info from a matrix",{
   mat <- matrix( c( 1,-2, 0
                   , 2, 0, 1
				   )
				, nrow=2
				, byrow=TRUE 
	#			, dimnames=list(c("a", "b"), c("x","y", "z"))
				)
   ei <- as.data.frame(as.editmatrix(mat))
   expect_equal(ei$edit, c("x1 == 2*x2", "2*x1 + x3 == 0"))
   
   mat <- matrix( c( 1,-2
                   , 2, 0
		 		       )
				    , nrow=2
				    , byrow=TRUE
                , dimnames=list(c("A", "B"), c("x","y"))
				)
   ei <- as.data.frame(as.editmatrix(mat))
   #expect_equal(ei$name, c("A","B"))
   expect_equal(ei$edit, c("x == 2*y", "2*x == 0"))
})

test_that("getb works",{
   cond <- c( "x + y > 2"
            , "y < 10"
            )
   E <- editmatrix(cond, FALSE)
   b <- getb(E)
  
   expect_equal(b, c(num1=2,num2=10))
})

test_that("getOps works",{
   cond <- c( "x + y > 2"
            , "y < 10"
            , "x + y == 2"
            , "y <= 10"
            , "y >= 10"
            )
   E <- editmatrix(cond, FALSE)
   ops <- getOps(E)
   expect_equivalent(ops, c(">","<","==","<=",">="))
})
