library(testthat)

context("Editmatrix")

test_that("editmatrix works correcly with character",{
   cond <- c( "x == y"
            , "z + w == y + x"
			   , "x + z == y + 2*w"
			   )
			
   mat <- editmatrix(cond)
   mat <- getA(mat)
   expect_equivalent(mat[1,], c(1,-1,0,0))
   expect_equivalent(mat[2,], c(-1,-1,1,1))
   expect_equivalent(mat[3,], c(1,-1,-2,1))
})

test_that("editmatrix works correcly with expression",{
  cond <- expression( x == y
             , z + w == y + x
             , x + z == y + 2*w
             )
  
  mat <- editmatrix(cond)
  mat <- getA(mat)
  expect_equivalent(mat[1,], c(1,-1,0,0))
  expect_equivalent(mat[2,], c(-1,-1,1,1))
  expect_equivalent(mat[3,], c(1,-1,-2,1))
})

test_that("editmatrix can simplify",{
   cond <- c( "2*x == x + y"
            , "z + 2 == y + 3"
            , "w == 3"
			   )
			
   E <- editmatrix(cond)
   mat <- getA(E)
   C <- getb(E)
   expect_equal(as.integer(mat[1,]), c(1,-1,0,0))
   expect_equal(as.integer(mat[2,]), c(0,-1,1,0))
   expect_equal(C[2], c(num2=1))
   
   expect_equal(mat[3,], c(x=0,y=0,z=0,w=1))
   expect_equal(C[3], c(num3=3))
})


test_that("editmatrix works correcly with data.frame",{
   
   edtinf.csv <- 
"name,edit
A,x == y
B,z + w == y + x
C,z == y + 2*w
"
   edtinf <- read.csv((con <- textConnection(edtinf.csv)))
   close(con)
			
   mat <- editmatrix(edtinf)
   A <- getA(mat)
   expect_equivalent(A[1,], c(1,-1,0,0))
   expect_equivalent(A[2,], c(-1,-1,1,1))
   expect_equivalent(A[3,], c(0,-1,-2,1))
})

test_that("editmatrix works with constants",{
   cond <- c( "x + y > 2"
            , "y < 10"
            )
   E <- editmatrix(cond)
   mat <- getA(E)
   expect_equal(as.integer(mat[1,]), c(-1,-1))
   expect_equal(as.integer(mat[2,]), c(0,1))
   expect_equal(as.integer(getb(E)), c(-2,10))
})

test_that("conditional statement parsing is not working..",{
   expect_error(editmatrix("if(x < 2) y > 4"))
})

test_that("editmatrix works with negative constants",{
   cond <- c( "x + y > -2"
            , "y < -10"
            )
   E <- editmatrix(cond)
   mat <- getA(E)
   expect_equal(as.integer(mat[1,]), c(-1,-1))
   expect_equal(as.integer(mat[2,]), c(0,1))
   expect_equal(as.integer(getb(E)), c(2,-10))
})

test_that("editmatrix works with negative coefficients",{
   cond <- c( "-2*x + y > 2"
            )
   E <- editmatrix(cond)
   mat <- getAb(E)
   expect_equivalent(mat[1,], c(2,-1,-2))
})

test_that("editmatrix works with coefficient after variable",{
  cond <- c( "x*-2 + y > 2"
             )
  E <- editmatrix(cond)
  mat <- getAb(E)
  expect_equivalent(mat[1,], c(2,-1,-2))
})

test_that("editmatrix fails with nonconstant coefficient",{
   cond <- c( "a*x == 2"
            )
   expect_error(editmatrix(cond))
})


test_that("is.editmatrix works",{
   mat <- editmatrix("x==y")
   expect_true(is.editmatrix(mat))
   expect_false(is.editmatrix(unclass(mat)))
})

test_that("as.editmatrix works",{
   A <- matrix( c( 1,-2, 0
                   , 2, 0, 1
				   )
				, nrow=2
				, byrow=TRUE
	#			, dimnames=list(c("a", "b"), c("x","y", "z"))
				)
   E <- as.editmatrix(A, b=c(0,1), ops=c("==","<"))
   ei <- as.data.frame(E)
   expect_equivalent(ei$edit, c("x1 == 2*x2", "2*x1 + x3 < 1"))
})

test_that("editmatrix normalize works",{
   cond <- c( "x > y"
            , "z + w >= y + x"
			   , "x + z < y + 2*w"
			   , "x + z == y + 2*w"
			   , "x + z <= y + 2*w"
			   )
			
   E <- editmatrix(editrules=cond, normalize=TRUE)
   mat <- getA(E)
   
   expect_equivalent(mat[1,], c(-1,1,0,0))
   expect_equivalent(mat[2,], c(1,1,-1,-1))
   expect_equivalent(mat[3,], c(1,-1,-2,1))
   expect_equivalent(mat[4,], c(1,-1,-2,1))
   expect_equivalent(mat[5,], c(1,-1,-2,1))
   ops <- getOps(E)
   expect_equivalent(ops, c("<", "<=", "<","==", "<="))
})


test_that("coercions work",{
    E <- editmatrix("x+y==z")
    expect_that(E, is_identical_to(editmatrix(as.data.frame(E)))) 
    expect_that(E, is_identical_to(as.editmatrix(A=getA(E), b=getb(E), ops=getOps(E))))    
    # edge case, testing as.character feature
    E <- editmatrix("x + 0.1*y==z")
    expect_that(E, is_identical_to(editmatrix(as.character(E)))) 
})

