context('Check filter_actors')

load('./tabari-demo.RData')

test_that("map actors", {
  dd <- map_actors(tabari.demo, list(osg=c('OSG'), gon.mor=c('GON', 'MOR')))
  expect_equal(nrow(dd), 66)

  line1 <- dd[1,]
  expect_equal(as.Date('1995-01-01'), line1$date)
  expect_equal("ARN", sources(line1))
  expect_equal("gon.mor", targets(line1))
  expect_equal("064", codes(line1))

  line3 <- dd[3,]
  expect_equal(as.Date('1995-01-02'), line3$date)
  expect_equal("gon.mor", sources(line3))
  expect_equal("osg", targets(line3))
  expect_equal("064", codes(line3))

  line60 <- dd[60,]
  expect_equal(as.Date('1995-01-14'), line60$date)
  expect_equal("gon.mor", sources(line60))
  expect_equal("ITH", targets(line60))
  expect_equal("066", codes(line60))
})



