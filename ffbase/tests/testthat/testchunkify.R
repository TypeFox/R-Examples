library(testthat)
library(ff)

context("chunkify")

test_that("Chunkify a function",{
	x <- 1:10 
	xf <- ff(x)
   
  sin.ff <- chunkify(sin)
   
	expect_identical( sin(x)
	                , sin.ff(xf)[]
				    )
})

test_that("Chunkexpr works with large expressions",{
  ep <- expression(symbol == "000001" & date >= as.Date("1993-03-05") & date <= as.Date("1993-04-13"))
  nm <- c("date", "symbol")
  e <- ffbase:::chunkexpr(ep, nm)
  expect_equal(length(e), 1)
})
