context("statcan")

##  Test the statcan algorithm
test <- read.csv("statcan.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
  expect_true(statcan(test$word[i]) == test$value[i])
test$test <- statcan(test$word)
for(i in 1:nrow(test))
  expect_true(test$test[i] == test$value[i])
