context("soundex")

##  Test the soundex algorithm
test <- read.csv("soundex.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
    expect_true(soundex(test$word[i]) == test$value[i])
test$test <- soundex(test$word)
for(i in 1:nrow(test))
    expect_true(test$test[i] == test$value[i])

##  Test the refined soundex algorithm
test <- read.csv("soundex-refined.csv", comment.char="#", stringsAsFactors=FALSE)
for(i in 1:nrow(test))
    expect_true(refinedSoundex(test$word[i]) == test$value[i])
test$test <- refinedSoundex(test$word)
for(i in 1:nrow(test))
    expect_true(test$test[i] == test$value[i])
