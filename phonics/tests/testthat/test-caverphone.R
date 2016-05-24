context("caverphone")

##  Test the Caverphone algorithm
test <- read.csv("caverphone.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
  expect_true(caverphone(test$word[i]) == test$value[i])
test$test <- caverphone(test$word)
for(i in 1:nrow(test))
  expect_true(test$test[i] == test$value[i])

##  Test the Caverphone 2 algorithm
test <- read.csv("caverphone-modified.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
  expect_true(caverphone(test$word[i], modified = TRUE) == test$value[i])
test$test <- caverphone(test$word, modified = TRUE)
for(i in 1:nrow(test))
  expect_true(test$test[i] == test$value[i])
