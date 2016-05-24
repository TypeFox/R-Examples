context("metaphone")

##  Test the metaphone algorithm.  Data comes from the Javascript
##  implementation.
test <- read.csv("metaphone.csv", comment.char="#", stringsAsFactors = FALSE)
for(i in 1:nrow(test))
    expect_true(metaphone(test$word[i]) == test$value[i])
test$test <- metaphone(test$word)
for(i in 1:nrow(test))
    expect_true(test$test[i] == test$value[i])
