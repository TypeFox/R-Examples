
today_pin <- substr(paste(paste(unlist(strsplit(as.character(Sys.Date()),split = "-")), collapse = ""),"0000",sep=""), 3, 12)
tomorrow_pin <- paste(paste(unlist(strsplit(as.character(Sys.Date()+1),split = "-")), collapse = ""),"0000",sep="")

context("as.pin")

test_that(desc="class: pin",{
  expect_is(as.pin("196408233234"), class = "pin")
})

test_that(desc="numeric: YYYYMMDDNNNC",{
  expect_equal(as.character(suppressMessages(as.pin(196408233234))), expected = "196408233234")
  expect_equal(as.character(suppressMessages(as.pin(200108230000))), expected = "200108230000")  
  expect_is(as.character(suppressMessages(as.pin(200108230000))), "character")  
  expect_is(as.character(suppressMessages(as.pin(pin = c(NA,198501169885)))), "character")
  expect_is(as.character(suppressMessages(as.pin(pin = as.numeric(NA)))), "character")
})

test_that(desc="numeric: YYMMDDNNNC",{
  expect_equal(as.character(suppressMessages(as.pin(6408233234))), expected = "196408233234")
  expect_equal(as.character(suppressMessages(as.pin(108230000))), expected = "200108230000")
  expect_equal(as.character(suppressMessages(as.pin(pin = c(8230000,108230000)))), expected = c("200008230000", "200108230000"))
  expect_is(as.character(suppressMessages(as.pin(c(NA,8501169885)))), "character")
})

test_that(desc="character: 'YYMMDDNNNC'",{
  expect_equal(as.character(suppressMessages(as.pin("6408233234"))), expected = "196408233234")  
  expect_equal(as.character(suppressMessages(as.pin("0008230000"))), expected = "200008230000")
  expect_equal(as.character(suppressMessages(as.pin(today_pin))), expected = paste("20",today_pin, sep=""))  
  expect_is(as.character(suppressMessages(as.pin(c(NA,"8501169885")))), "character")
})

test_that(desc="factor: 'YYMMDDNNNC'",{
  expect_equal(as.character(suppressMessages(as.pin(as.factor("6408233234")))), expected = "196408233234")  
  expect_equal(as.character(suppressMessages(as.pin(as.factor("0008230000")))), expected = "200008230000")
})

test_that(desc="character: 'YYYYMMDDNNNC'",{
  expect_equal(as.character(suppressMessages(as.pin("196408233234"))), expected = "196408233234")  
  expect_is(as.character(suppressMessages(as.pin(c(NA,"198501169885")))), "character")
})

test_that(desc="different formats",{
  expect_equal(as.character(suppressMessages(as.pin(c("196408233234", "640823-3234", "19640823-3234", "6408233234")))), 
                          expected = rep("196408233234", 4))
  expect_equal(as.character(suppressMessages(as.pin(c(196408233234, 6408233234)))), 
               expected = rep("196408233234", 2))
})

test_that(desc="character: 'YYMMDD-NNNC'",{
  expect_equal(as.character(as.pin("640823-3234")), expected = "196408233234")  
  expect_equal(as.character(as.pin("000823-0000")), expected = "200008230000")
  expect_equal(as.character(as.pin("000823+0000")), expected = "190008230000")
})

test_that(desc="error expected",{
  suppressWarnings(expect_equal(as.character(as.pin(tomorrow_pin)), as.character(NA)))
  suppressWarnings(expect_equal(as.character(as.pin(pin = "AA6408233234")), as.character(NA)))
  suppressWarnings(expect_equal(as.character(as.pin("196418233234")), as.character(NA)))
  suppressWarnings(expect_equal(as.character(as.pin("196408333234")), as.character(NA)))
  expect_warning(as.pin(tomorrow_pin))
  expect_warning(as.pin("AA6408233234"))
  expect_warning(as.pin("196418233234"))
  expect_warning(as.pin("196408333234"))
  
  test_pin <- c("196408233234", tomorrow_pin, "AA6408323234", "19640823323", "1964083332349", "196408333234", "19640823-334", "19640823")
  test_pin_res <- c(TRUE, rep(FALSE, 7))
  suppressWarnings(expect_equal(!is.na(as.pin(test_pin)), test_pin_res))
  
  non_relevant_class <- lm(1:10~rep(1:5,2))
  expect_error(as.pin(non_relevant_class))
  expect_error(as.pin(c(TRUE,FALSE)))
})


test_that("as.pin.pin", {
  suppressWarnings(expect_equal(as.pin(as.pin("test_pin")), as.pin("test_pin")))
})

test_that("as.pin.logical", {
  expect_is(as.pin(NA), "pin")
  expect_error(as.pin(TRUE))
})


test_pins <- c("18920822-2298", "18920822-2299", "19920419-1923")
test_that("Recycling rules", {
  expect_is(data.frame(as.pin(test_pins), 1:9), "data.frame")
  expect_equal(nrow(data.frame(as.pin(test_pins), 1:9)), 9)
  expect_equal(data.frame(as.pin(test_pins), 1:9)[1:3, 1], data.frame(as.pin(test_pins), 1:9)[4:6, 1])
  expect_equal(data.frame(as.pin(test_pins), 1:9)[1:3, 1], data.frame(as.pin(test_pins), 1:9)[7:9, 1])
})


semi_pins <- c("550504333A", "19280118123X", "850504111T", "850504111 ", "19280118123 ")
test_that("deceased 1947 - 1967", {
  suppressWarnings(expect_is(as.pin(semi_pins), "pin"))
  expect_warning(as.pin(semi_pins[3]), "Erroneous pin")
  expect_warning(as.pin(semi_pins[4]), "Erroneous pin")
  expect_message(as.pin(semi_pins[1]), "less than 100 years old and people with birth year")
  expect_message(as.pin(semi_pins[2]), "Assumption: People with birth year before 1967 and character")
  expect_message(as.pin(semi_pins[5]), "Assumption: People with birth year before 1967 and character")
})


test_that("Expect message only when YYMMDDNNNC format is used", {
  num_to_check <- c("202100-6255","121212-1212","19121212-1212","121212+1212", 121212121212, NA, Inf, TRUE, F, "foo", 123, 456L)
  expect_silent(suppressWarnings(as.pin(num_to_check)))
})
