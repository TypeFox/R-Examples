test_that("Ping API", {    
  expect_that(PingAPI(), equals("Success"))
})