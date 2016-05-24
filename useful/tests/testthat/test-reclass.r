context('Checking that reclass sets the class accordingly')

# create data.frames and lists to test with
df1 <- data.frame(A=1:10, B=1:10)
list1 <- list(A=1:3, B=2)
class(list1) <- c('Special', 'list')
df2 <- df1
list2 <- list1

# reclass them using both methods
reclass(df1) <- 'NewClass'
df2 <- reclass(df2, 'NewClass')
reclass(list1) <- 'NewClass'
list2 <- reclass(list2, 'NewClass')

test_that('reclass returns the appropriate number of classes', {
    expect_equal(length(class(df1)), 2)
    expect_equal(length(class(df2)), 2)
    expect_equal(length(class(list1)), 3)
    expect_equal(length(class(list2)), 3)
})

test_that('reclass returns the correct classes', {
    expect_equal(class(df1), c('NewClass', 'data.frame'))
    expect_equal(class(df2), c('NewClass', 'data.frame'))
    expect_equal(class(list1), c('NewClass', 'Special', 'list'))
    expect_equal(class(list2), c('NewClass', 'Special', 'list'))
})