context("String manipulation")

test_getBetweenRight <- function(corrChars,
                                 left = '(',
                                 right = ')',
                                 testString = "([1234&5])",
                                 insertSt = NULL){

  mapply(function(rt, outChars){
    parse <- getBetween(testString, left, rt, insertSt)
    expect_equal(nchar(parse), outChars,
      info = paste0(left, ' with ', rt, '->', parse))
  }, right, corrChars)
}

test_getBetweenLeft <- function(corrChars,
                                left = '(',
                                right = ')',
                                testString = "([1234&5])",
                                insertSt = NULL){
  mapply(function(lt, outChars){
    parse <- getBetween(testString, lt, right, insertSt)
    expect_equal(nchar(parse), outChars,
      info = paste0(lt, ' with ', right, '->', parse))
  }, left, corrChars)
}

test_that("getBetween does croping right with 1 symbol", {
  test_getBetweenLeft(c(8, 7, 2), left = c("(", "[", "&"))
  test_getBetweenRight(c(7, 5), right = c("]", "&"))
})

test_that("getBetween does croping right with 2 symbols", {
  test_getBetweenLeft(c(7, 6, 1), left = c("(\\[", "[1", "&5"))
  test_getBetweenRight(c(7, 5), right = c("]\\)", "&5"))
})

test_that("getBetween does insert right with 1 symbol", {
  test_getBetweenLeft(c(6, 7, 12), left = c("(", "[", "&"), insertSt = "four")
  test_getBetweenRight(c(6, 7), right = c(")", "]"), insertSt = "four")
})

test_that("getBetween does insert right with 2 symbols", {
  test_getBetweenLeft(c(7, 8, 13), left = c("(\\[", "[1", "&5"), insertSt = "four")
  test_getBetweenRight(c(7, 9), right = c("]\\)", "&5"), insertSt = "four")
})

test_that("getBetween does one character default",{
  test_getBetweenLeft(c(1), left = "")
  test_getBetweenRight(1, right = "")
  expect_match(getBetween("asdf(%1)", '%', ''),
    '1')
  expect_equal(match(getBetween("asd(%1)sd", '%', '', 'c'),"asd(%c)sd"),
    1)
})


test_that("getBetween does inclusion right", {
  res <- getBetween("sdf[1234]jkly", left = "[", right = "]", shInclude = TRUE)
  expect_true(grepl("\\[1234\\]", res))

  res <- getBetween("sdf[1234]jkly",
    left = "[", right = "]", insertSt <- "hello", shInclude = TRUE)
  expect_true(grepl("sdfhellojkly", res))
  
  res <- getBetween("sdf[www1234]wwwjkly",
  	left = "[www", right = "]www", insertSt <- "hello", shInclude = TRUE)
  expect_true(grepl("sdfhellojkly", res))

})

test_that("getBetween picks up nothing if one of the positions aren't found",{
	res <- getBetween("dfg}sdfg", ".", "}")
	expect_true(grepl(res, ""))
	
	
})
test_that("getBetween inserts nothing if one of the positions aren't found",{
	res <- getBetween("dfg}sdfg", ".", "}", insertChar = "123123")
	expect_equal(match(res, "dfg}sdfg"),1)
	
})

test_that("getBetween is vectorized",{
	res <- getBetween(rep("{dfg}sdfg", 2), "{", "}")
	expect_true(all(grepl("dfg",res)))
	
})
