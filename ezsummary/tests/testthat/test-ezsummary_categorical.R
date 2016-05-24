context("Summary_Categorical")
#####

df <- mtcars %>% group_by(am) %>% select(gear)
df2 <- mtcars %>% group_by(am) %>% select(gear, carb)
df3 <- mtcars %>% select(gear)

test_that("ezsummary_categorical can work without grouping data", {
  expect_equal(names(ezsummary_categorical(df3)), expected = c("variable", "count", "p"))
  expect_equal(names(ezsummary_categorical(df3, n=T)), expected = c("variable", "N", "count", "p"))
  expect_equal(ezsummary_categorical(df3)[2], expected = structure(list(count = c(15L, 12L,5L)), class = c("tbl_df", "data.frame"), row.names = 1:3))
  expect_equal(ezsummary_categorical(df3, n=T)[2:3], expected = structure(list(N=c(32L, 32L, 32L), count = c(15L, 12L,5L)), class = c("tbl_df", "data.frame"), row.names = 1:3))
})

test_that("ezsummary_categorical can rename the variables for data with grouping info", {
  expect_equal(names(ezsummary_categorical(df2))[2], expected = c("variable"))
  expect_equal(names(ezsummary_categorical(df2, n=T))[2:3], expected = c("variable", "N"))
  expect_equivalent(as.character(ezsummary_categorical(df2)[2, 2]), expected = c("gear_4"))
  expect_equivalent(as.character(ezsummary_categorical(df2, n=T)[2, 2:3]), expected = c("gear_4", "19"))
})

test_that("ezsummary_categorical can export n.group and n.var", {
  expect_equal(attributes(ezsummary_categorical(df))$n.group, expected = 1)
  expect_equal(attributes(ezsummary_categorical(df))$n.var, expected = 1)
  expect_equal(attributes(ezsummary_categorical(df2))$n.group, expected = 1)
  expect_equal(attributes(ezsummary_categorical(df2))$n.var, expected = 2)
})
