context('Checking that upper and lower case are detected properly.')

# create text to check
toCheck <- c('BIG', 'little', 'Mixed', 'BIG WITH SPACE', 'little with space', 'Mixed With Space', 
             'UPPER 17', 'UPPER17', 'UPP17ER', 'lower 17', 'lower17', 'low17er',
             'Mixed 17', 'Mixed17', 'Mix17ed',
             '19')

test_that('The functions return logicals', {
    expect_is(find.case(toCheck, 'upper'), 'logical')
    expect_is(find.case(toCheck, 'lower'), 'logical')
    expect_is(find.case(toCheck, 'mixed'), 'logical')
    expect_is(find.case(toCheck, 'numeric'), 'logical')
    
    expect_is(upper.case(toCheck), 'logical')
    expect_is(lower.case(toCheck), 'logical')
    expect_is(mixed.case(toCheck), 'logical')
})

test_that('The correct number of elements are returned', {
    expect_equal(length(find.case(toCheck, 'upper')), length(toCheck))
    expect_equal(length(find.case(toCheck, 'lower')), length(toCheck))
    expect_equal(length(find.case(toCheck, 'mixed')), length(toCheck))
    expect_equal(length(find.case(toCheck, 'numeric')), length(toCheck))
    
    
    expect_equal(length(upper.case(toCheck)), length(toCheck))
    expect_equal(length(lower.case(toCheck)), length(toCheck))
    expect_equal(length(mixed.case(toCheck)), length(toCheck))
    
    expect_equal(length(find.case(toCheck, 'upper')), length(upper.case(toCheck)))
    expect_equal(length(find.case(toCheck, 'lower')), length(lower.case(toCheck)))
    expect_equal(length(find.case(toCheck, 'mixed')), length(mixed.case(toCheck)))
})

test_that('find.case returns the same results as upper.case and lower.case', {
    expect_identical(find.case(toCheck, 'upper'), upper.case(toCheck))
    expect_identical(find.case(toCheck, 'lower'), lower.case(toCheck))
    expect_identical(find.case(toCheck, 'mixed'), mixed.case(toCheck))
})

test_that('find.case, upper.case and lower.case work on simple data', {
    expect_identical(find.case(toCheck, 'upper'), c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
                                                    TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
                                                    FALSE, FALSE, FALSE,
                                                    TRUE))
    expect_identical(upper.case(toCheck), c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
                                            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
                                            FALSE, FALSE, FALSE,
                                            TRUE))
    
    expect_identical(find.case(toCheck, 'lower'), c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
                                                    FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
                                                    FALSE, FALSE, FALSE,
                                                    TRUE))
    expect_identical(lower.case(toCheck), c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
                                            FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
                                            FALSE, FALSE, FALSE,
                                            TRUE))
    
    expect_identical(find.case(toCheck, 'mixed'), c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
                                                    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                    TRUE, TRUE, TRUE,
                                                    FALSE))
    expect_identical(mixed.case(toCheck), c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
                                            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                            TRUE, TRUE, TRUE,
                                            FALSE))
    
    expect_identical(find.case(toCheck, 'numeric'), c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                    FALSE, FALSE, FALSE,
                                                    TRUE))
})