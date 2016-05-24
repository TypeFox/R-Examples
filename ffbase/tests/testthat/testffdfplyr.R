library(testthat)
library(ff)

context("ffdfplyr")

test_that("ffdfplyr works",{
	data(iris)
	ffiris <- as.ffdf(iris)
	testFUN <- function(x){
		lowestbypetalwidth <- x[order(x$Petal.Width, decreasing=TRUE), ]
 		lowestbypetalwidth <- lowestbypetalwidth[!duplicated(lowestbypetalwidth[, c("Species","Petal.Width")]), ]
 		lowestbypetalwidth$group <- factor(x= "lowest", levels = c("lowest","highest"))
 		highestbypetalwidth <- x[order(x$Petal.Width, decreasing=FALSE), ]
 		highestbypetalwidth <- highestbypetalwidth[!duplicated(highestbypetalwidth[, c("Species","Petal.Width")]), ]
 		highestbypetalwidth$group <- factor(x= "highest", levels = c("lowest","highest"))
 		rbind(lowestbypetalwidth, highestbypetalwidth)
	}
	test.ff <- ffdfdply(x = ffiris, split = ffiris$Species, FUN = function(x) testFUN(x), BATCHBYTES = 5000, trace=TRUE)
	test.ram <- split(iris, iris$Species)
	test.ram <- lapply(test.ram, FUN=function(x) testFUN(x))
	test.ram <- do.call(rbind, test.ram)

	expect_true(!sum(!apply(test.ff[,], MARGIN=1, FUN=function(x) paste(x, collapse=",")) %in% apply(test.ram, MARGIN=1, FUN=function(x) paste(x, collapse=","))))
	expect_true(!sum(!apply(test.ram, MARGIN=1, FUN=function(x) paste(x, collapse=",")) %in% apply(test.ff[,], MARGIN=1, FUN=function(x) paste(x, collapse=","))))
})



