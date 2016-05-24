test_that("test.is_cas_number.a_character_vector.returns_true_when_string_contains_a_cas_number", 
  {
    x <- c(
      water = "7732-18-5", d_glucose = "50-99-7", 
      l_glucose = "921-60-8", no_hyphens = "7732185", 
      two_check_digits = "7732-18-55", bad_check_digit = "7732-18-4"
    )
    expected <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
    expect_equal(strip_attributes(actual <- is_cas_number(x)), expected)
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "bad format", "bad checkdigit"), c(3, 2, 1)))
    )
  })

test_that("test.is_credit_card_number.invalid_card_numbers.returns_false_for_all", 
  {
    x <- c(
      "4012888888881882", 
      "4111 1111 1111 11111", "5655 5555 5555 4443", 
      "51051 051 0510 5100", "3782 822463 1005"
    )
    expected <- rep.int(FALSE, length(x))
    expect_equal(strip_attributes(actual <- is_credit_card_number(x)), expected)
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("bad checkdigit", "bad format"), c(1, 4)))
    )
  })

test_that("test.is_credit_card_number.valid_card_numbers.returns_true_for_all", 
  {
    x <- c("4111 1111 1111 1111", "4012888888881881", "5555 5555 5555 4444", 
      "5105 1051 0510 5100", "3782 822463 10005", "3714 496353 98431", 
      "3787 344936 71000", "3056 930902 5904", "3852 000002 3237", "6011 1111 1111 1117", 
      "6011 0009 9013 9424", "3530 1113 3330 0000", "3566 0020 2036 0505")
    expected <- rep.int(TRUE, length(x))
    names(expected) <- x
    expect_equal(is_credit_card_number(x), expected)
  })

test_that("test.is_email_address.a_character_vector_rfc5322_match.returns_true_when_string_contains_an_email_address", 
{
  x <- c("foo@bar.com", "foo@@bar.com", "@bar.com", "foo@bar", "foo@bar.comma", 
    "foo!@bar.com", NA)
  expected <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, NA)
  expect_equal(
    strip_attributes(actual <- is_email_address(x, method = "rfc5322")), 
    expected
  )
  expect_equal(names(actual), unname(x))
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "bad format", "", "missing"), c(1, 3, 2, 1)))
  )
})

test_that("test.is_email_address.a_character_vector_simple_match.returns_true_when_string_contains_an_email_address", 
{
  x <- c("foo@bar.com", "foo@@bar.com", "@bar.com", "foo@bar", "foo@bar.comma", 
         "foo!@bar.com", NA)
  expected <- rep.int(c(TRUE, FALSE, NA), c(1, 5, 1))
  expect_equal(
    strip_attributes(actual <- is_email_address(x)), 
    expected
  )
  expect_equal(names(actual), unname(x))
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "bad format", "missing"), c(1, 5, 1)))
  )
})

test_that("test.is_hex_color.a_character_vector.returns_true_when_string_contains_a_hex_colour", 
{
  x <- c("#123456", "#789aBc", "#dEF000", "123456", "g12345", "#12 34 56", 
    "#12345", "#1234567", NA)
  expected <- rep.int(c(TRUE, FALSE, NA), c(3, 5, 1))
  expect_equal(
    strip_attributes(actual <- is_hex_color(x)), 
    expected
  )
  expect_equal(names(actual), unname(x))
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "bad format", "missing"), c(3, 5, 1)))
  )
})

test_that("test.is_ip_address.a_character_vector.returns_true_when_string_contains_an_ip_address", 
{
  x <- c(localhost = "localhost", valid_address = "255.0.255.0", out_of_range = "1.2.3.256", 
    five_blocks = "1.2.3.4.5", non_numeric = "1.2.3.Z", missing_block = "1.2.3.NA", NA)
  expected <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_ip_address(x)), 
    expected
  )
  expect_equal(names(actual), unname(x))
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "big numbers", "bad format", "missing"), c(2, 1, 3, 1)))
  )
})

test_that("test.is_isbn10_code.a_character_vector_type_10.returns_true_when_string_contains_an_isbn10_code", 
{
  x <- c(hyphens = "0-387-98503-4", spaces = "0 387 98503 4", just_numbers = "0387985034", 
    too_long = "00-387-98503-4", too_short = "0-387-9850-4", non_numeric = "Z-387-98503-4", 
    invalid_check_digit = "0-387-98503-5", missing = NA)
  expected <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_isbn10_code(x)), 
    expected
  )
  expect_equal(names(actual), unname(x))
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "bad length", "bad format", "bad check digit", "missing"), c(3, 2, 1, 1, 1)))
  )
})

test_that("test.is_isbn13_code.a_character_vector_type_13.returns_true_when_string_contains_an_isbn13_code", 
  {
    x <- c(hyphens = "978-0-387-98503-9", spaces = "978 0 387 98503 9", 
      just_numbers = "9780387985039", too_long = "9978-0-387-98503-9", 
      too_short = "978-0-387-9850-9", non_numeric = "Z78-0-387-9850-9", 
      invalid_check_digit = "978-0-387-98503-8", missing = NA)
    expected <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, NA)
    expect_equal(
      strip_attributes(actual <- is_isbn13_code(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "bad length", "bad format", "bad check digit", "missing"), c(3, 2, 1, 1, 1)))
    )
  })
