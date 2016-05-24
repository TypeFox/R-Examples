context('encapsulation')

tst <- lrequire(encapsulation)

# test1.1 -- test to see if variables exposed in the lrequire-d file exist in the current namespace

test_that('lrequired files encapsulate non-exposed variables', {
  expect_false(exists('will.not.be.exposed'))
})

test_that('exposed variable is set and known', {
  expect_equal(tst$hello.world, 'hello world!')
})

test_that('exposed function is callable and functions correctly', {
  expect_equal(tst$sum(1, 2), 3)
})
