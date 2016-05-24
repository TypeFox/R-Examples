test_that('is.algebra', {
    for (a in names(.algebras)) {
        expect_true(is.algebra(algebra(a)))
    }
})
