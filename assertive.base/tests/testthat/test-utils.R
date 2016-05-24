test_that(
  "to_names with a character vector returns unquoted character names",
  {
    x <- c("abc", "", "NA", NA)
    expected <- x
    actual <- to_names(x)
    expect_equal(actual, expected)
  }
)

test_that(
  "to_names with a (double) numeric vector returns 17 sig fig 'g' formatting",
  {
    x <- with(.Machine, c(1 + double.eps, double.xmax, double.xmin, NaN, -Inf, NA))
    expected <- ifelse(is.na(x), NA_real_, sprintf("%.17g", x))
    actual <- to_names(x)
    expect_equal(actual, expected)
  }
)

test_that(
  "to_names with a complex vector returns returns 17 sig fig 'g' formatting for each component",
  {
    x <-  with(.Machine, c(1 + double.eps, double.xmax, double.xmin, NaN, -Inf, NA)) * (1 + 1i)
    expected <- ifelse(is.na(x), NA_complex_, sprintf("%.17g+%.17gi", Re(x), Im(x)))
    actual <- to_names(x)
    expect_equal(actual, expected)
  }
)

test_that(
  "to_names with a list input returns deparsed object",
  {
    x <- list(NA, 1 + .Machine$double.eps, "abc", list(Inf, list("def", NA)))
    expected <- as.character(x)
    actual <- to_names(x)
    expect_equal(actual, expected)
  }
)



test_that("test.coerce_to.numeric_vector_to_data_frame.returns_data_frame", 
  {
    x <- 1:5
    expected <- data.frame(x = x)
    expect_equal(suppressWarnings(coerce_to(x, "data.frame")), expected)
    expect_warning(coerce_to(x, "data.frame"))
  })

test_that(
  "dont_stop with multiple errors and warnings successfully runs",
  {
    expected <- list(
      'stop("If you don\'t stop;")' = simpleError("If you don't stop;"),
      'warning("Someone\'s gonna find yo\' ass dead (this is a warning)")' = simpleWarning("Someone's gonna find yo' ass dead (this is a warning)"),
      'warning("Someone\'s gonna poison your food (this is a warning)")' = simpleWarning("Someone's gonna poison your food (this is a warning)"),
      'stop("Don\'t stop, no no, you\'ll be sorry")' = simpleError("Don\'t stop, no no, you\'ll be sorry"),
      'stop("Don\'t stop, thinking about tomorrow")' = simpleError("Don't stop, thinking about tomorrow"),
      'stop("Don\'t stop, it\'ll soon be here")' = simpleError("Don't stop, it'll soon be here")
    )
    actual <- dont_stop(
      {
        # With apologies to Lil' Kim
        stop("If you don't stop;")
        warning("Someone's gonna find yo' ass dead (this is a warning)")
        warning("Someone's gonna poison your food (this is a warning)")
        stop("Don't stop, no no, you'll be sorry")
        
        # Bonus errors for David, Jenny and other Fleetwood Mac fans
        stop("Don't stop, thinking about tomorrow")
        stop("Don't stop, it'll soon be here")
      }
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "dont_stop works with objects that don't deparse to a single string",
  {
    expected <- list("function() {}" = function() {})
    actual <- dont_stop(
      # deparse returns a character vector 
      function() {} 
    )
    # don't test for identicality due to function environment
    expect_equal(actual, expected)
  }
)

test_that(
  "get_name_in_parent works when object exists outside of function",
  {
    outside <- 1
    f <- function(inside) 
    {
      get_name_in_parent(inside)
    }
    expected <- "outside"
    actual <- f(outside)
    expect_identical(actual, expected)
  }
)

test_that(
  "get_name_in_parent works when object doesn't exist outside of function",
  {
    f <- function(inside) 
    {
      get_name_in_parent(inside)
    }
    expected <- "1"
    actual <- f(1)
    expect_identical(actual, expected)
  }
)

test_that(
  "get_name_in_parent works with percent inside function call, escape_percent = TRUE",
  {
    f <- function(inside) 
    {
      get_name_in_parent(inside)
    }
    expected <- "1 %%>%% exp"
    actual <- f(1 %>% exp)
    expect_identical(actual, expected)
  }
)

test_that(
  "get_name_in_parent works with percent inside function call, escape_percent = FALSE",
  {
    f <- function(inside) 
    {
      get_name_in_parent(inside, FALSE)
    }
    expected <- "1 %>% exp"
    actual <- f(1 %>% exp)
    expect_identical(actual, expected)
  }
)

test_that(
  "merge_dots_with_list takes duplicates from ...",
  {
    expected <- list(x = 1, y = "b", z = TRUE)
    expect_warning(
      actual <- merge_dots_with_list(x = 1, y = "b", l = list(y = "c", z = TRUE)),
      "Duplicated arguments: y"
    )
    expect_equal(actual, expected)
  }
)

test_that(
  "merge_dots_with_list throws an error with unnamed arguments and allow_unnamed_elements = FALSE",
  {
    expect_error(
      merge_dots_with_list(x = 1, "b", l = list(y = "c")),
      "There are unnamed elements in x or y, but allow_unnamed_elements = FALSE"
    )
  }
)

test_that(
  "merge_dots_with_list works with unnamed arguments and allow_unnamed_elements = TRUE",
  {
    expected <- list(x = 1, y = "c", "b")
    actual <- merge_dots_with_list(
      x = 1, "b", l = list(y = "c"), 
      allow_unnamed_elements = TRUE
    )
    expect_equal(actual, expected)
  }
)

test_that(
  "merge_dots_with_list works with no list argument",
  {
    expected <- list(x = 1, y = "b")
    actual <- merge_dots_with_list(x = 1, y = "b")
    expect_equal(actual, expected)
  }
)

test_that(
  "test.parenthesise.character_input.returns_parenthesised_input",  
  {
    x <- "foo"
    types <- eval(formals(parenthesise)$type)
    actual <- vapply(
      types,
      function(type) parenthesise(x, type),
      character(1),
      USE.NAMES = FALSE
    )
    expected <- c(
      "(foo)", "[foo]", "{foo}", "<foo>", "\u3008foo\u3009", 
      "- foo -", "\u2013 foo \u2013", "\u2014foo\u2014", ", foo, "
    )
    expect_identical(actual, expected)
  }
)

test_that("test.use_first.a_list_double_indexing.returns_contents_of_first_element", 
  {
    x <- as.list(letters)
    expected <- "a"
    expect_identical(expected, suppressWarnings(use_first(x)))
    expect_warning(use_first(x))
  })

test_that("test.use_first.a_list_single_indexing.returns_first_element", {
  x <- as.list(letters)
  expected <- list("a")
  expect_identical(expected, suppressWarnings(use_first(x, "[")))
  expect_warning(use_first(x, "["))
})

test_that("test.use_first.a_scalar.returns_x", {
  x <- "a"
  expected <- x
  expect_identical(expected, use_first(x))
})

test_that("test.use_first.a_vector_double_indexing.returns_first_element", 
  {
    x <- letters
    expected <- "a"
    expect_identical(expected, suppressWarnings(use_first(x)))
    expect_warning(use_first(x))
  })

test_that("test.use_first.a_vector_single_indexing.returns_first_element", 
  {
    x <- letters
    expected <- "a"
    expect_identical(expected, suppressWarnings(use_first(x, "[")))
    expect_warning(use_first(x, "["))
  })

test_that("test.use_first.empty.throws_error", {
  x <- NULL
  expect_error(use_first(x))
}) 
