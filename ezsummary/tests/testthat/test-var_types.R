context("var_type")
#####

df1 <- mtcars

test_that("var_types will only react to data.frame(including data.table)", {
  expect_error(var_types(NULL, "c"), "Please supply a data.frame/data.table as the value of tbl")
  expect_error(var_types("ABC", "c"), "Please supply a data.frame/data.table as the value of tbl")
  expect_error(var_types(1:10, "c"), "Please supply a data.frame/data.table as the value of tbl")
  expect_error(var_types(list(1,2), "c"), "Please supply a data.frame/data.table as the value of tbl")
  expect_error(var_types(list(1:10,"a"), "c"), "Please supply a data.frame/data.table as the value of tbl")
})

test_that("var_types will throw error when the length of 'types' doesn't match up", {
  expect_error(var_types(mtcars, "qqqqqqc"), "The length of the string you entered doesn't match the number of variables")
  expect_error(var_types(group_by(mtcars, am), "qcqqqqqccc"), "The length of the string you entered doesn't match the number of variables")
})

test_that("var_types will throw error when the 'types' string contains characters other than 'q' and 'c'", {
  expect_error(var_types(mtcars, "qcqqqqqcccw"), 'Unrecognizable character\\(s\\) detected!! Please review your input and use "q" and "c" to denote quantitative and categorical variables')
})

test_that("var_types can add the desired contents to the attribute of the data.frame", {
  expect_equal(attributes(var_types(mtcars, "qcqqqqqcccc"))$var_types, expected = c("q", "c", "q", "q", "q", "q", "q", "c", "c", "c", "c"))
})
