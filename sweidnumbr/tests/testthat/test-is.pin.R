
tomorrow_pin <- paste(paste(unlist(strsplit(as.character(Sys.Date()+1),split = "-")), collapse = ""),"0000",sep="")
test_pin <- c("196408233234")

context("is.pin")

test_that(desc="is.pin",{
  expect_is(is.pin(pin = test_pin), "logical")
  expect_equal(is.pin(pin = test_pin), expected = FALSE)
})
