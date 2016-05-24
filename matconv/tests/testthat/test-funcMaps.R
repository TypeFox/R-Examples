context("Function Maps")


test_that("Basic argument switching", {
  dict <- "matsort:sort, 2, 1"
  example <- c("asdf", "hjkl")
  map <- makeFuncMaps(dict)
  result <- map$matsort$argMap[[1]](example)$rargs
  expect_match(result, paste(rev(example), collapse = ', '))
})

test_that("Literal numbers in output", {
  dict <- "matsort:sort, 2, 1L"
  example <- c("asdf", "hjkl")
  map <- makeFuncMaps(dict)
  result <- map$matsort$argMap[[1]](example)$rargs
  expect_match(result, paste(example[2], 1, sep = ', '))
})

test_that("Literal arg inserts", {
  dict <- "matsort:sort, %1 * %2, asdf%1"
  example <- c("asdf", "hjkl")
  map <- makeFuncMaps(dict)
  result <- map$matsort$argMap[[1]](example)$rargs
  expect_equal(result,
    paste(
      paste0("sort(", example[1], ' * ', example[2]),
      paste0("asdf", example[1]),
      sep = ', '))
})

test_that("Flag switching integrates", {
	dict <- c(
		"matsort--if 2:sort, 2, 1",
		"matsort--if 3:sort, 3")
	example <- c(
		"thing <- matsort(asdf, hjkl)",
		"thing <- matsort(asdf, hjkl, omg)")
	map <- makeFuncMaps(dict)
	result <- convFunctionsCalls(example, map)
	expect_equal(result[1],
		"thing <- sort(hjkl, asdf)")
	expect_equal(result[2],
		"thing <- sort(omg)")

})

test_that("Multiple outputs work", {
	dict <- "matsort:sort, 2, 1 --out mean std"
	example <- "[myMean myStd] <- matsort(asdf, hjkl)"
	map <- makeFuncMaps(dict)
	result <- convFunctionsCalls(example, map)
	expect_equal(result[1],
		"lout <- sort(hjkl, asdf); myMean <- lout$mean; myStd <- lout$std")

})

test_that("Can parse using space separated args", {
	dict <- "matsort:sort, 2, 1 --space-sep"
	example <- "matsort asdf hjkl"
	map <- makeFuncMaps(dict)
	result <- convFunctionsCalls(example, map)
	expect_equal(result[1],
		"sort(hjkl, asdf)")

})

test_that("Can convert functions within the line", {
	dict <- c(
		"matsort--if 2:sort, 2, 1",
		"matsort--if 3:sort, 3")
	example <- c(
		"thing <- 6 * 9 * lifeOffset(matsort(asdf, hjkl))")
	map <- makeFuncMaps(dict)
	result <- convFunctionsCalls(example, map)
	expect_equal(result[1],
		"thing <- 6 * 9 * lifeOffset(sort(hjkl, asdf))")
	
})

test_that("Can convert functions with strings in it", {
	dict <- c(
		"matsort--if 2:sort, 2, 1",
		"matsort--if 3:sort, 3")
	example <- c(
		"thing <- 6 * 9 * lifeOffset(matsort('Thing', hjkl))")
	map <- makeFuncMaps(dict)
	result <- convFunctionsCalls(example, map)
	expect_equal(result[1],
		"thing <- 6 * 9 * lifeOffset(sort(hjkl, 'Thing'))")
	
})
