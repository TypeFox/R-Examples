context("format.pin")

x <- structure(c("191212121212", "201212121212", "191212121212", "201212121212", 
            "199106252523", "189611070798", "190011121298", "200809050523", 
            "189101252598", "201401232323"), class = c("AsIs", "pin", "character"))

suppressMessages({
  
  x <- as.pin(fake_pins$pin)
  test_that("format.pin works", {
    expect_is(format(x), "character")
    expect_equal(format_pin(x[1], "%N"), "1212")
    expect_equal(format_pin(x[1], "%P"), "(19) 12-12-12 - 1212")
    expect_equal(format_pin(x[1], "%Y-%m-%d-%N"), "1912-12-12-1212")
  })
  
})