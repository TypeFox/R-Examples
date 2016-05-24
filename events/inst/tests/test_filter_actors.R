context('Check filter_actors')

load('./tabari-demo.RData')

test_that("test both", {
  dd <- filter_actors(tabari.demo, spotter('GON','OSG'))
  expect_equal(nrow(dd), 7)

  dd2 <- filter_actors(tabari.demo, spotter('GON','OSG'), which='both')
  expect_identical(dd, dd2)

  line1 <- dd[1,]
  expect_equal(as.Date('1995-01-02'), line1$date)
  expect_equal("GON", sources(line1))
  expect_equal("OSG", targets(line1))
  expect_equal("064", codes(line1))

  line2 <- dd[7,]
  expect_equal(as.Date('1995-01-05'), line2$date)
  expect_equal("GON", sources(line2))
  expect_equal("OSG", targets(line2))
  expect_equal("141", codes(line2))
    
})

test_that("test source", {
  dd <- filter_actors(tabari.demo, spotter('GON','OSG'), which='source')
  expect_equal(nrow(dd), 25)

  line1 <- dd[1,]
  expect_equal(as.Date('1995-01-02'), line1$date)
  expect_equal("GON", sources(line1))
  expect_equal("OSG", targets(line1))
  expect_equal("064", codes(line1))

  line2 <- dd[25,]
  expect_equal(as.Date('1995-01-17'), line2$date)
  expect_equal("GON", sources(line2))
  expect_equal("BRE", targets(line2))
  expect_equal("033", codes(line2))

})

test_that("test target", {
  dd <- filter_actors(tabari.demo, spotter('GON','OSG'), which='target')
  expect_equal(nrow(dd), 26)

  line1 <- dd[1,]
  expect_equal(as.Date('1995-01-01'), line1$date)
  expect_equal("ARN", sources(line1))
  expect_equal("GON", targets(line1))
  expect_equal("064", codes(line1))

  line2 <- dd[26,]
  expect_equal(as.Date('1995-01-17'), line2$date)
  expect_equal("BRE", sources(line2))
  expect_equal("GON", targets(line2))
  expect_equal("032", codes(line2))
})



