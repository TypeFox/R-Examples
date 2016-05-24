context('Checking that classdf returns the right type, the correct number of elements')

test_that('classdf returns lists and vectors as appropriate', {
    expect_is(classdf(CO2), 'list')
    expect_is(classdf(iris), 'character')
    expect_is(classdf(mtcars), 'character')
})

test_that('classdf returns the correct number of elements', {
    expect_equal(length(classdf(CO2)), length(CO2))
    expect_equal(length(classdf(iris)), length(iris))
    expect_equal(length(classdf(mtcars)), length(mtcars))
})

test_that('classdf returns a named object', {
    expect_named(classdf(CO2))
    expect_named(classdf(iris))
    expect_named(classdf(mtcars))
})

test_that('classdf returns an error if data is not a data.frame', {
    expect_error(classdf(list(A=1:3, B=2)))
})