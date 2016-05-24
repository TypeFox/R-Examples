context("Data Conversion")


test_that("Can turn into matrix: Digits",{
	
	exm <- "1, 5,6;2,7,2"
	matRep <- matrixify(exm)
	expect_true(grepl("nrow = 2", matRep))
	expect_true(grepl("ncol = 3", matRep))
	expect_true(grepl("c\\(1, 5, 6, 2, 7, 2\\)", matRep))
})

test_that("Can turn into matrix: Numerics",{
	
	exm <- "1, 5,6;2,7,2.234"
	matRep <- matrixify(exm)
	expect_true(grepl("nrow = 2", matRep))
	expect_true(grepl("ncol = 3", matRep))
	expect_true(grepl("c\\(1, 5, 6, 2, 7, 2.234\\)", matRep))
})

test_that("Integration: data maps", {
	dataMap <- makeDataMap("[", "]", "matrix")
	
	test <- "scrRes <- [23,2, 3.2; 7, 6, 8];"
	res <- dataMap(test)
	eval(parse(text = res))
	
	expect_true(is.matrix(scrRes))
	expect_equal(dim(scrRes)[1], 2)
	expect_equal(dim(scrRes)[2], 3)
	
	dataVecMap <- makeDataMap("{", "}", "vector")
	testCell <- "celRes <- {23,2, 3.2; 7, 6, 8};"
	res <- dataVecMap(testCell)
	eval(parse(text = res))
	
	expect_true(is.vector(celRes))
	expect_equal(length(celRes), 6)
	expect_null(dim(celRes))
	
})

test_that("Data maps doesn't pick up non-instantiations",{
	dataMap <- makeDataMap("{", "}", "matrix")
	negTest <- "stmTable{rr,4} <- 'Feed'"
	res <- dataMap(negTest)
	
	expect_true(!is.na(match(res, negTest)))
	
	
})

test_that("Missing numbers can be inputted as data",{
	dataMap <- makeDataMap("{", "}", "matrix")
	NATest <- "stmTable <- {NaN}"
	res <- dataMap(NATest)
	eval(parse(text = res))
	expect_true(is.nan(stmTable[1]))
})

test_that("Convert string data",{
	dataMap <- makeDataMap(rClass = "vector", matClass = "string")
	
	test <- "csvFile <- [beg num2str(ispc()) '.csv']"
	res <- dataMap(test)
	
	expect_equal(1, match("csvFile <- c(paste0(beg, num2str(ispc()), '.csv'))", res))
})


test_that("Convert string data with escaped quotes",{
	dataMap <- makeDataMap(rClass = "vector", matClass = "string")
	
	test <- "comd <- [realWd '/xls2csv.vbs \"' fullpath '\" ' int2str(numSheetsParsed)]"
	res <- dataMap(test)
	
	expect_true(grepl("paste", res))
	expect_true(grepl("xls2csv", res))
	
	test2 <- "comd <- ['ssconvert -S \"' fullpath '\" \"' beg '%n.csv\"']"
	res2 <- dataMap(test2)
	expect_equal(1, match(res2, 
		"comd <- c(paste0('ssconvert -S \"', fullpath, '\" \"', beg, '%n.csv\"'))"))
})

test_that("Does not convert data as a string concat",{
	dataMap <- makeDataMap(rClass = "vector", matClass = "string")
	
	test <- "comd <- [23,2, 3.2; 7, 6, 8]"
	res <- dataMap(test)
	
	expect_true(match(test, res) == 1)
})

test_that("Does not convert string cat as matrix", {
	dataMap <- makeDataMap(matClass = "matrix", rClass = "matrix")
	
	test <- "comd <- [realWd '/xls2csv.vbs \"' fullpath '\" ' int2str(numSheetsParsed)]"
	res <- dataMap(test)
	
	expect_true(match(test, res) == 1)
	
	
	
})

test_that("R lists can be selected", {
	dataMap <- makeDataMap(matClass = "cell", rClass = "list")
	test <- "comd <- {23,2, 4; 'string', 6, 8}"
	res <- dataMap(test)
	
	expect_equal("comd <- list(list(23, 2, 4), list('string', 6, 8))", res)
	
})
