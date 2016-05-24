context("Slicing Data")


test_that("Can slice with one instance",{
	matLine <- "asdfasdf <- dat{1,2}"
	mp <- makeSliceMap("{", "}", "list")
	res <- mp(matLine)
	expect_true(!is.na(match(res, "asdfasdf <- dat[[1,2]]")))
	
	
})

test_that("Can slice with Structures",{
	matLine <- "asdfasdf <- dat.temp(1,2)"
	mp <- makeSliceMap(matClass = "structure", rClass = "list")
	res <- mp(matLine)
	expect_true(!is.na(match(res, "asdfasdf <- dat[['temp']](1,2)")))
	
	
})

test_that("Can slice with mult instance",{
	matLine <- "asdfasdf <- dat{3,4}{1,2}"
	mp <- makeSliceMap("{", "}", "list")
	res <- mp(matLine)
	expect_true(!is.na(match(res, "asdfasdf <- dat[[3,4]][[1,2]]")))
	
	
})

test_that("Can slice with mult struct instances",{
	matLine <- "asdfasdf <- dat.temp.K"
	mp <- makeSliceMap(matClass = "structure", rClass = "list")
	res <- mp(matLine)
	expect_true(!is.na(match(res, "asdfasdf <- dat[['temp']][['K']]")))
	
	
})

test_that("Negative case doesn't do anything",{
	matLine <- "asdfas <- dat"
	mp <- makeSliceMap(matClass = "structure", rClass = "list")
	res <- mp(matLine)
	expect_true(!is.na(match(res, "asdfas <- dat")))
})

test_that("Doesn't pick up cases in strings",{
	matLine <- "csvFile <- [beg num2str(sh - isunix()) 'h.csv']"
	mp <- makeSliceMap(matClass = "structure", rClass = "list")
	res <- mp(matLine)
	expect_true(!is.na(match(res, matLine)))
	
})
