context("Hetu (personal identification number in Finland")

test_that("hetu works correctly", {
  expect_equal(hetu(c("010101-0101", "111111-111C"), extract = "gender"), c("Female", "Male"))
  expect_equal(hetu(c("111111-111C"), extract = "gender"), "Male")
  expect_equal(as.character(hetu("010101-0101")$hetu), "010101-0101")
  expect_equal(as.character(hetu(c("010101-0101"))$gender), "Female")
  expect_equal(hetu("010101-0101")$personal.number, 10)
  expect_equal(as.character(hetu("010101-0101")$checksum), "1")
  expect_equal(as.character(hetu("010101-0101")$date), "1901-01-01")
  expect_equal(hetu("010101-0101")$day, 1)
  expect_equal(hetu("010101-0101")$month, 1)
  expect_equal(hetu("010101-0101")$year, 1901)
  expect_equal(as.character(hetu("010101-0101")$century.char), "-")

  expect_true(valid_hetu("010101-0101"))
  expect_false(valid_hetu("010101-010A"))
})