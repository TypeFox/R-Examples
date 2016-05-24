library(testthat)

context("append")

test_that("Appending ff vector works",{
   x <- c(1,3)
   y <- c(2,3)
   fx <- ff(x)
   fy <- ff(y)
   
   fz <- ffappend(NULL, x)
   expect_equal(length(fz), length(x))
   fz <- ffappend(fz, y)
   expect_equal(length(fz), length(x)+length(y))
   expect_identical(fz[], c(x,y))
   
   fz <- ffappend(NULL, fx)
   expect_equal(length(fz), length(fx))
   fz <- ffappend(fz, fy)
   expect_equal(length(fz), length(fx)+length(fy))
   expect_identical(fz[], c(x,y))
})

test_that("Concatenating ff vector works",{
   x <- ff(c(1,3))
   y <- ff(c(2,3))
   z <- c(x)
   expect_equal(length(z), length(x))
   z <- c(x, y)
   expect_equal(length(z), length(x)+length(y))
   expect_identical(z[], c(x[],y[]))
   z <- c(x, y, x)
   expect_equal(length(z), length(x)+length(y)+length(x))

   a <- ff(factor(LETTERS[1:5]))
   b <- ff(factor(LETTERS[5:7]))
   expect_identical(c(as.character(a[]), as.character(b[])), as.character(ffappend(a, b)[]))

})

test_that("Appending dataframe to ffdf works",{
   dat <- data.frame(x=1:3, y=3:1, z=as.factor(c("a","b","c")))
	fdat <- NULL
	fdat <-ffdfappend(fdat, dat)
	expect_identical(dat,fdat[,])
	fdat <- ffdfappend(fdat,dat)
	expect_equal(nrow(fdat),2*nrow(dat))
})

test_that("Appending ffdf to ffdf works",{
   dat <- as.ffdf(data.frame(x=1:3, y=3:1, z=as.factor(c("a","b","c"))))
   fdat <- NULL
   fdat <-ffdfappend(fdat, dat)
   expect_identical(dat[,],fdat[,])
   fdat <- ffdfappend(fdat,dat)
   expect_equal(nrow(fdat),2*nrow(dat))
})

test_that("Coercing an ff vector works",{
   x <- 1:100
   y <- rnorm(100)
   fx <- ff(x)
   fy <- ff(y)
   
   fz <- ffappend(fx, fy, by = 2)
   expect_identical(fz[], c(x,y))
   
   x <- rep(NA, 100)
   y <- rnorm(100)
   fx <- ff(x)
   fy <- ff(y)
   
   fz <- ffappend(fx, fy, by = 2)
   expect_identical(fz[], c(x,y))
   
   x <- rep(TRUE, 100)
   y <- rep(NA, 100)
   fx <- ff(x, vmode = "boolean")
   fy <- ff(y)
   
   fz <- ffappend(fx, fy, by = 2)
   expect_identical(fz[], c(x,y))
 })
 
test_that("Coercing an ffdf works",{
  dat <- data.frame(x=rnorm(5), y=as.integer(1:5), z=as.factor(c("a","b","c","d","e")))
  x <- ffdf(x=ff(as.logical(rep(NA, 100))), y=ff(as.numeric(rnorm(100))), z=ff(as.factor(rep("Z", 100))))

	fdat <- ffdfappend(x=x, dat=dat, by = 2)
	expect_identical(rbind(x[,], dat),fdat[,])
})

test_that("Appending two ffdfs works",{
  a <- data.frame(a=1:10, b=factor(rep(c("A", "B"), 5)))
  b <- as.ffdf(a)
  c <- as.ffdf(a)
  expect_equivalent(ffdfappend(b, c)[,], rbind(a,a))
})
