
today_pin <- substr(paste(paste(unlist(strsplit(as.character(Sys.Date()),split = "-")), collapse = ""),"0000",sep=""), 3, 12)
tomorrow_pin <- substr(paste(paste(unlist(strsplit(as.character(Sys.Date()+1),split = "-")), collapse = ""),"0000",sep=""), 3, 12)

context("pin_century")

test_that(desc="century",{
  expect_equal(pin_century(pin = today_pin), expected = 20)
  expect_equal(pin_century(pin = tomorrow_pin), expected = 19)
  expect_is(pin_century(pin = today_pin), "numeric")
})

