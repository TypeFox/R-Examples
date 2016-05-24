context("mra")

##  This is structured a bit differently from the other tests in order
##  to test both the encoder and the comparison.

##  Test the MRA encoding algorithm
test <- read.csv("mra-encode.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
    expect_true(mra_encode(test$word[i]) == test$value[i])
test$test <- mra_encode(test$word)
for(i in 1:nrow(test))
    expect_true(test$test[i] == test$value[i])

##  Test the MRA compare algorithm
test <- read.csv("mra-compare.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
    expect_true(mra_compare(mra_encode(test$word1[i]), mra_encode(test$word2[i])) == test$value[i])
test$test <- mra_compare(mra_encode(test$word1), mra_encode(test$word2))
for(i in 1:nrow(test))
    expect_true(test$test[i] == test$value[i])
