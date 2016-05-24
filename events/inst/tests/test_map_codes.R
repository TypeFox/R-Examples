context('Check filter_actors')

load('./tabari-demo.RData')

test_that("map codes", {
  dd <- map_codes(tabari.demo, list(sixtyfour='064', 
     fortyone.thirtytwo=c('041', '032')))
  line1 <- dd[1,]
  expect_equal(as.Date('1995-01-01'), line1$date)
  expect_equal("ARN", sources(line1))
  expect_equal("GON", targets(line1))
  expect_equal("sixtyfour", codes(line1))

  line63 <- dd[63,]
  expect_equal(as.Date('1995-01-17'), line63$date)
  expect_equal("BRE", sources(line63))
  expect_equal("GON", targets(line63))
  expect_equal("fortyone.thirtytwo", codes(line63))
})
