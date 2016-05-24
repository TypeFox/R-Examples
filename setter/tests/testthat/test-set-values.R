context("test set/copy class")

test_that(
  "set_class works",
  {
    expected <- structure(1:12, class = "foo")
    expect_warning(
      actual <- 1:12 %>% set_class(factor("foo")),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_class works",
  {
    expected <- as.list(1:12)
    actual <- 1:12 %>% copy_class(list())
    expect_identical(actual, expected)
  }
)

context("test set/copy mode")

test_that(
  "set_mode works",
  {
    expected <- 3L
    expect_warning(
      actual <- pi %>% set_mode(factor("integer")),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_mode works",
  {
    expected <- (1:12) + 0
    actual <- 1:12 %>% copy_mode(pi)
    expect_identical(actual, expected)
  }
)

context("test set/copy storage_mode")

test_that(
  "set_storage_mode works",
  {
    expected <- 3L
    expect_warning(
      actual <- pi %>% set_storage_mode(factor("integer")),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_storage_mode works",
  {
    expected <- (1:12) + 0
    actual <- 1:12 %>% copy_storage_mode(pi)
    expect_identical(actual, expected)
  }
)

context("test set/copy dim")

test_that(
  "set_dim works",
  {
    expected <- matrix(1:12, 3)
    expect_warning(
      actual <- 1:12 %>% set_dim(c(3, 4)),
      "Coercing value to class .integer."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_dim works",
  {
    expected <- matrix(1:12, 3)
    actual <- 1:12 %>% copy_dim(matrix(raw(12), 3))
    expect_identical(actual, expected)
  }
)

context("test set/copy length")

test_that(
  "set_length works",
  {
    expected <- c(1:12, NA)
    expect_warning(
      actual <- 1:12 %>% set_length(13),
      "Coercing value to class .integer."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_length works",
  {
    expected <- c(1:12, NA)
    actual <- 1:12 %>% copy_length(1:13)
    expect_identical(actual, expected)
  }
)

context("test set/copy names")

test_that(
  "set_names works",
  {
    expected <- structure(1:12, names = month.name)
    expect_warning(
      actual <- 1:12 %>% set_names(factor(month.name)),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_names works",
  {
    expected <- structure(1:12, names = month.name)
    actual <- 1:12 %>% copy_names(structure(raw(12), names = month.name))
    expect_identical(actual, expected)
  }
)

context("test set/copy colnames")

test_that(
  "set_colnames works",
  {
    expected <- data.frame(January = 1:5, February = pi ^ (1:5))
    expect_warning(
      actual <- data.frame(x = 1:5, y = pi ^ (1:5)) %>%
        set_colnames(factor(month.name[1:2])),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_colnames works",
  {
    expected <- data.frame(January = 1:5, February = pi ^ (1:5))
    actual <- data.frame(x = 1:5, y = pi ^ (1:5)) %>%
      copy_colnames(data.frame(January = 1, February = 2))
    expect_identical(actual, expected)
  }
)

context("test set/copy rownames")

test_that(
  "set_rownames works with a matrix",
  {
    expected <- matrix(1:12, 3, dimnames = list(month.name[1:3], NULL))
    expect_warning(
      actual <- matrix(1:12, 3) %>% set_rownames(factor(month.name[1:3])),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "set_rownames works with a data.frame & factor names",
  {
    expected <- data.frame(x = 1:5, y = pi ^ (1:5), row.names = LETTERS[1:5])
    expect_warning(
      actual <- data.frame(x = 1:5, y = pi ^ (1:5)) %>%
        set_rownames(factor(LETTERS[1:5])),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "set_rownames works with a list converted to a data.frame & missing names",
  {
    expected <- data.frame(x = 1:5, y = pi ^ (1:5))
    actual <- list(x = 1:5, y = pi ^ (1:5)) %>%
      set_class("data.frame") %>%
      set_rownames()
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_rownames works",
  {
    expected <- data.frame(x = 1:5, row.names = LETTERS[1:5])
    actual <- data.frame(x = 1:5) %>%
      copy_rownames(data.frame(y = pi ^ (1:5), row.names = LETTERS[1:5]))
    expect_identical(actual, expected)
  }
)

context("test set/copy dimnames")

test_that(
  "set_dimnames works",
  {
    expected <- matrix(1:12, 3, dimnames = list(letters[1:3], LETTERS[1:4]))
    actual <- matrix(1:12, 3) %>%
      set_dimnames(list(letters[1:3], LETTERS[1:4]))
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_dimnames works",
  {
    expected <- matrix(1:12, 3, dimnames = list(letters[1:3], LETTERS[1:4]))
    actual <- matrix(1:12, 3) %>%
      copy_dimnames(data.frame(A = 1:3, B = 4:6, C = 7:9, D = 10:12, row.names = letters[1:3]))
    expect_identical(actual, expected)
  }
)

context("test set/copy levels")

test_that(
  "set_levels works",
  {
    expected <- factor(letters[1:5], letters[5:1])
    actual <- 5:1 %>%
      set_levels(letters[5:1]) %>%
      set_class("factor")
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_levels works",
  {
    expected <- factor(letters[5:1], letters[5:1])
    actual <- factor(letters[1:5]) %>%
      copy_levels(factor(sample(letters[1:5], 10, replace = TRUE), letters[5:1]))
    expect_identical(actual, expected)
  }
)

context("test set/copy comment")

test_that(
  "set_comment works",
  {
    expected <- structure(list(), comment = as.character(1:10))
    expect_warning(
      actual <- list() %>% set_comment(1:10),
      "Coercing value to class .character."
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_comment works",
  {
    expected <- structure(TRUE, comment = LETTERS)
    actual <- TRUE %>%
      copy_comment(structure(list(), comment = LETTERS))
    expect_identical(actual, expected)
  }
)

# context("test set/copy formals")
#
# test_that(
#   "set_formals works",
#   {
#     skip("TODO")
#     expected <- function(z, ...) {sin(pi)}
#     actual <- function(x, y = 2, ...) {sin(pi)} %>%
#         set_formals(alist(z = , ... = ))
#     expect_identical(actual, expected)
#   }
# )
#
# test_that(
#   "copy_formals works",
#   {
#     skip("TODO")
#     expected <- structure(TRUE, comment = LETTERS)
#     actual <- TRUE %>%
#       copy_formals(structure(list(), comment = LETTERS))
#     expect_identical(actual, expected)
#   }
# )
#
# context("test set/copy body")
#
# test_that(
#   "set_body works",
#   {
#     skip("TODO")
#   }
# )
#
# test_that(
#   "copy_body works",
#   {
#     skip("TODO")
#   }
# )

context("test set/copy attributes")

test_that(
  "set_attributes works",
  {
    expected <- structure(1:12, foo = pi, bar = exp(1), baz = 99)
    expect_warning(
      actual <- 1:12 %>%
        set_attributes(foo = pi, bar = exp(1), .dots = list(bar = exp(2), baz = 99)),
      "Duplicated arguments: bar"
    )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_attributes works",
  {
    expected <- structure(1:12, foo = pi, baz = 99)
    actual <- 1:12 %>%
      copy_attributes(
        structure(FALSE, foo = pi, bar = exp(1), baz = 99),
        c("foo", "baz")
      )
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_all_attributes works",
  {
    expected <- structure(1:12, foo = pi, bar = exp(1), baz = 99)
    actual <- 1:12 %>%
      copy_all_attributes(structure(FALSE, foo = pi, bar = exp(1), baz = 99))
    expect_identical(actual, expected)
  }
)

test_that(
  "copy_most_attributes works",
  {
    # dim and dimnames shouldn't be copied
    expected <- structure(list(1, 2:3), foo = pi)
    actual <- list(1, 2:3) %>%
      copy_most_attributes(structure(matrix(FALSE, dimnames = list("a", "A")), foo = pi))
    expect_identical(actual, expected)
  }
)



