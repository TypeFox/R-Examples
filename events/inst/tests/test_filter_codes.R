context('Check filter_codes')

load('./tabari-demo.RData')

test_that("test codes", {
  dd <- filter_codes(tabari.demo, spotter('064', '041'))
  expect_equal(nrow(dd), 25)

  line1 <- dd[1,]
  expect_equal(as.Date('1995-01-01'), line1$date)
  expect_equal("ARN", sources(line1))
  expect_equal("GON", targets(line1))
  expect_equal("064", codes(line1))

  line2 <- dd[25,]
  expect_equal(as.Date('1995-01-04'), line2$date)
  expect_equal("ORC", sources(line2))
  expect_equal("OSG", targets(line2))
  expect_equal("041", codes(line2))
})



