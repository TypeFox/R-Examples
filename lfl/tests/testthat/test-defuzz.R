test_that('defuzz mom', {
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.9, 0.2, 0), 1:10, type='mom'),
                 7)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.5, 0.9, 0.2, 0), 1:10, type='mom'),
                 7)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.8, 0.2, 0), 1:10, type='mom'),
                 6.5)
})

test_that('defuzz fom', {
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.9, 0.2, 0), 1:10, type='fom'),
                 6)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.5, 0.9, 0.2, 0), 1:10, type='fom'),
                 6)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.8, 0.2, 0), 1:10, type='fom'),
                 6)
})

test_that('defuzz lom', {
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.9, 0.2, 0), 1:10, type='lom'),
                 8)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.5, 0.9, 0.2, 0), 1:10, type='lom'),
                 8)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.8, 0.2, 0), 1:10, type='lom'),
                 7)
})

test_that('defuzz dee', {
    # as mom
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.9, 0.2, 0), 1:10, type='dee'),
                 7)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.5, 0.9, 0.2, 0), 1:10, type='dee'),
                 7)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.8, 0.2, 0), 1:10, type='dee'),
                 6.5)

    # as fom
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.9, 0.9, 0.9), 1:10, type='dee'),
                 6)
    expect_equal(defuzz(c(0, 0, 0, 0.1, 0.3, 0.9, 0.9, 0.9, 0.9, 1), 1:10, type='dee'),
                 10)

    # as lom
    expect_equal(defuzz(c(0.9, 0.9, 0.9, 0.9, 0.9, 0.7, 0.5, 0, 0, 0), 1:10, type='dee'),
                 5)
    expect_equal(defuzz(c(1, 0.9, 0.9, 0.9, 0.9, 0.7, 0.5, 0, 0, 0), 1:10, type='dee'),
                 1)
})
