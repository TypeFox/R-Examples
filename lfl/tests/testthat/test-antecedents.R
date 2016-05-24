test_that('antecedents', {
    expect_equal(antecedents(list(c('a', 'b', 'c'),
                                  c('d'),
                                  c('a', 'e'))),
                 list(c('b', 'c'),
                      character(0),
                      c('e')))
})

test_that('antecedents (single rule)', {
    expect_equal(antecedents(list(c('a', 'b', 'c'))),
                 list(c('b', 'c')))
})

