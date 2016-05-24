context("Easy syntax")

test_that("End within slicing gets replaced", {
	linesMat <- c("for(ca in commas(2:end))")
	res <- convSymbols(linesMat)
	expect_true(!is.na(match(res, "for(ca in commas(2:length(commas)))")))
	
	linesMat <- c("for(ca in commas(2:end)){")
	res <- convSymbols(linesMat)
	expect_true(!is.na(match(res, "for(ca in commas(2:length(commas))){")))
	
})

test_that("End to get the last element gets replaced", {
	
	linesMat <- c("thing{end}")
	res <- convSymbols(linesMat)
	expect_true(!is.na(match(res, "thing{length(thing)}")))
})
