context("cologne")

##  Test the Cologne algorithm
test <- read.csv("cologne.csv", comment.char="#", stringsAsFactors = FALSE, colClasses = c("character", "character"), encoding = "UTF-8")
for(i in 1:nrow(test))
  expect_true(cologne(test$word[i]) == test$value[i])
test$test <- cologne(test$word)
for(i in 1:nrow(test))
  expect_true(test$test[i] == test$value[i])
