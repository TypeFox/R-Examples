test_that('consequents', {
    expect_equal(consequents(list(c('a', 'b', 'c'),
                                  c('d'),
                                  c('a', 'e'))),
                 list('a', 'd', 'a'))
})
