context("rogerroot")

##  Test the Roger Root name encoding
test <- read.csv("rogerroot.csv", comment.char="#", stringsAsFactors=FALSE, colClasses=c("character", "character"))
for(i in 1:nrow(test))
    expect_true(rogerroot(test$word[i]) == test$value[i])
test$test <- rogerroot(test$word)
for(i in 1:nrow(test))
    expect_true(test$test[i] == test$value[i])
