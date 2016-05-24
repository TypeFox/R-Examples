context('Checking that subspecials replaces need special characters with escaped characters.')

# Create fake data with missing values
needSubs <- c("Hello", "(parens)", "Excited! Mark")

# run command a bunch of times
out1 <- subOut(c("Hello", "(parens)", "Excited! Mark"))
out2 <- subOut(c("Hello", "(parens)", "Excited! Mark"), specialChars=c("!", "("))
specials1 <- subSpecials(c("Hello", "(parens)", "Excited! Mark"))
specials2 <- subSpecials(c("Hello", "(parens)", "Excited! Mark"), specialChars=c("!", "("))
specials3 <- subSpecials(c("Hello", "(parens)", "Excited! Mark"), c("This is a period. And this is an asterisk *"), specialChars=c("!", "("))
specials4 <- subSpecials(c("Hello", "(parens)", "Excited! Mark"), c("This is a period. And this is an asterisk *"), specialChars=c("!", "(", "*"))

# set some expectations
eo1 <- c("Hello", "\\(parens\\)", "Excited\\! Mark")
eo2 <- c("Hello", "\\(parens)", "Excited\\! Mark")
es1 <- list(c("Hello", "\\(parens\\)", "Excited\\! Mark"))
es2 <- list(c("Hello", "\\(parens)", "Excited\\! Mark"))
es3 <- list(c("Hello", "\\(parens)", "Excited\\! Mark"), "This is a period. And this is an asterisk *")
es4 <- list(c("Hello", "\\(parens)", "Excited\\! Mark"), "This is a period. And this is an asterisk \\*")

test_that('subOut returns the correct number of elements and type', {
    expect_is(out1, 'character')
    expect_is(out2, 'character')
    
    expect_equal(length(out1), 3)
    expect_equal(length(out2), 3)
})

test_that('subSpecials returns the correct number of elements and type', {
    expect_is(specials1, 'list')
    expect_is(specials2, 'list')
    expect_is(specials3, 'list')
    expect_is(specials4, 'list')

    expect_equal(length(specials1), 1)
    expect_equal(length(specials2), 1)
    expect_equal(length(specials3), 2)
    expect_equal(length(specials4), 2)
    
    expect_equal(length(specials1[[1]]), 3)
    expect_equal(length(specials2[[1]]), 3)
    expect_equal(length(specials3[[1]]), 3)
    expect_equal(length(specials3[[2]]), 1)
    expect_equal(length(specials4[[1]]), 3)
    expect_equal(length(specials4[[2]]), 1)
})

test_that('subOut operates on special characters', {
    expect_identical(out1, eo1)
    expect_identical(out2, eo2)
})

test_that('subSpecials operates on special characters', {
    expect_identical(specials1, es1)
    expect_identical(specials2, es2)
    expect_identical(specials3, es3)
    expect_identical(specials4, es4)
})