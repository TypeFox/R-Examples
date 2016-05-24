
context("pin_birthplace")

today_pin <- paste(paste(unlist(strsplit(as.character(Sys.Date()),split = "-")), collapse = ""),"0000",sep="")
pin_test <- c("0000000019876", "187001019876","196408233234", "196408833234", today_pin, "196408830000")
pin_test_res <-
  c(NA, "Extra number and immigrants (immigrated after 1946)",
    "Gotlands l\u00E4n", NA,
    "Born after 31 december 1989",
    NA)

test_that(desc="birthplace",{
  suppressWarnings(expect_equal(as.character(pin_birthplace(pin = pin_test)), expected = pin_test_res))
  suppressWarnings(expect_is(pin_birthplace(pin = pin_test), "factor"))
})

test_that(desc="Handle NA, interim and coordn in pin_birthplace",{
  suppressWarnings(expect_true(is.na(pin_birthplace(pin = as.pin(c("hejbaberiba","198501169885")))[1])))
  suppressWarnings(expect_true(all(is.na(pin_birthplace(pin = c("19000625P816","190006859816"))))))
})
