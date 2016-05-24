context("nysiis")

##  Test the NYSIIS algorithm
test <- read.csv("nysiis.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
  expect_true(nysiis(test$word[i]) == test$value[i])
test$test <- nysiis(test$word)
for(i in 1:nrow(test))
  expect_true(test$test[i] == test$value[i])

##  Test the modified NYSIIS algorithm
test <- read.csv("nysiis-modified.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
  expect_true(nysiis(test$word[i], modified = TRUE) == test$value[i])
test$test <- nysiis(test$word, modified = TRUE)
for(i in 1:nrow(test))
  expect_true(test$test[i] == test$value[i])
