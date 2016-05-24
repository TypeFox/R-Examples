
context("Adjustment methods")

library(editrules)
E <- editmatrix(expression(
   x + y == 1,
   x > 0,
   y >= 0)
)

test_that("adjust.editarray works dense",{

   expect_equal(adjust(E,c(x=0,y=0))$x, c(x=0.5,y=0.5))
   # no var labels
   expect_equal(adjust(E,c(0,0))$x, c(0.5,0.5))
   # more var in x than in E
   expect_equal(adjust(E,c(x=0,y=0,z=10))$x, c(x=0.5,y=0.5,z=10))
})

test_that("adjust.editarray works sparse",{

	expect_equal(
		adjust(E,c(x=0,y=0),method='dense')$x, 
		adjust(E,c(x=0,y=0),method='sparse')$x
	) 
	expect_equal(
		adjust(E,c(0,0),method='dense')$x, 
		adjust(E,c(0,0),method='sparse')$x
	) 
	expect_equal(
		adjust(E,c(x=0,y=0,z=0),method='dense')$x, 
		adjust(E,c(x=0,y=0,z=0),method='sparse')$x
	) 

})


test_that("adjust.matrix",{
   A <- getA(E)
   b <- getb(E)
   ops <- getOps(E)

   I <- order(ops,decreasing=TRUE)
   A <- A[I,,drop=FALSE]
   b <- b[I]
   neq <- sum(ops == "==")

   expect_equivalent(adjust(A,b,c(0,0),neq=neq)$x,c(0.5,0.5))
   expect_error( adjust(A,b,c(0,NA),neq=neq) )
   expect_error( adjust(A,b,c(0,0,0),neq=neq) )
})

test_that("adjust notes exceeding maxiter",{
   e <- editmatrix(expression(x < -1, x > 0))
   expect_true(grep("maximum",adjust(e,0)$status)==1)
})




