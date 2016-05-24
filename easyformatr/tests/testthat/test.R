test_that("filter_info", {
  expect_equal(
    filter_info(type == "time" & component == "base"),
    structure(c("|code |base         |", "|:----|:------------|",
                "|p    |am_pm        |", "|C    |century      |", "|F    |date         |",
                "|c    |datetime     |", "|d    |day          |", "|w    |day_of_week  |",
                "|j    |day_of_year  |", "|H    |hour         |", "|M    |minute       |",
                "|m    |month        |", "|n    |newline      |", "|S    |second       |",
                "|R    |time         |", "|s    |timestamp    |", "|Z    |timezone     |",
                "|V    |week_of_year |", "|t    |tab          |", "|Y    |year         |"
    ), format = "markdown", class = "knitr_kable") )
  expect_equal(
    filter_info(type == "number" & component == "base"),
    structure(c("|code |base       |", "|:----|:----------|", "|i    |integer    |",
                "|o    |octal      |", "|x    |hex        |", "|f    |double     |",
                "|e    |scientific |", "|g    |auto       |", "|a    |binary     |"
    ), format = "markdown", class = "knitr_kable") )
  expect_equal(
    filter_info(type == "time" & component == "mutant" & base == "second"),
    structure(c("|code |decimal | digits|", "|:----|:-------|------:|",
                "|0S   |TRUE    |     NA|", "|0S0  |NA      |      0|", "|0S1  |NA      |      1|",
                "|0S2  |NA      |      2|", "|0S3  |NA      |      3|", "|0S4  |NA      |      4|",
                "|0S5  |NA      |      5|", "|0S6  |NA      |      6|"), format = "markdown", class = "knitr_kable") )
  expect_equal(
    filter_info(type == "number" & component == "flag"),
    structure(c("|code |flag           |", "|:----|:--------------|",
                "|-    |left_justify   |", "|+    |always_sign    |", "|     |prefix_space   |",
                "|0    |zero_pad       |", "|#    |hex_prefix     |", "|#    |always_decimal |",
                "|#    |remove_zeros   |"), format = "markdown", class = "knitr_kable") )
  expect_equal(
    filter_info(type == "number" & component == "option"),
    structure(c("|option         |description                               |",
                "|:--------------|:-----------------------------------------|",
                "|before_decimal |number of digits before the decimal place |",
                "|after_decimal  |number of digits after the decimal place  |",
                "|use_input      |sprintf input argument number             |"
    ), format = "markdown", class = "knitr_kable") )
  expect_error(
    filter_info(component == "not a component"),
    "No information matches those conditions" )

  expect_equal(
    easy_format(easy_format("We are the 99%") )[1],
    "We are the 99%%")
})

test_that("easy_format", {
  expect_equal(
    easy_format(year, month, day, integer, octal, double)[1],
    "%Y%m%d%.i%.o%.f")
  expect_equal(
    easy_format(second %>% decimal)[1],
    "%0S")
  expect_equal(
    easy_format(list(integer,
                     double) %>%
                  always_decimal)[1],
    "%#.i%#.f")
  expect_equal(
    easy_format(list(month,
                     list(day,
                          minute) ) %>%
                  roman)[1],
    "%0m%0d%0M")
  expect_equal(
    easy_format(year, month, day, integer, octal, double)[1],
    "%Y%m%d%.i%.o%.f")
  expect_equal(
    easy_format(second %>% roman),
    easy_format(second) )
  expect_equal(
    easy_format(second %>% decimal),
    easy_format(second %>% decimal(TRUE) ) )
  expect_equal(
    easy_format(second %>%
                  decimal %>%
                  decimal(NA) ),
    easy_format(second) )
  expect_equal(
    easy_format(double %>% before_decimal),
    easy_format(double) )
  expect_equal(
    easy_format(double %>% before_decimal(3) %>% before_decimal),
    easy_format(double) )
  expect_equal(
    easy_format("We are the 99%")[1],
    "We are the 99%%")
  expect_equal(
    easy_format(roman("I am Spartacus"))[1],
    "I am Spartacus")
  expect_equal(
    easy_format("We", "are", "the 99%", sep = " ")[1],
    "We are the 99%%")
  expect_equal(
    easy_format(list(double %>%
                       zero_pad,
                     "1%",
                     double %>%
                       left_justify) %>%
                  use_input(value = 1) %>%
                  always_sign %>%
                  before_decimal(value = 3) %>%
                  after_decimal(value = 0),
                sep = " and ")[1],
    "%1$+03.0f and 1%% and %1$-+3.0f")
  expect_equal(
    easy_format(list(year %>% religious,
                     "/",
                     month %>% name) %>%
                  short,
                "/",
                list(day,
                     tab,
                     hour %>% twelve,
                     ":",
                     minute) %>%
                  roman,
                ":",
                second %>% digits(value = 3),

                " ",
                am_pm,
                " ",
                timezone)[1],
    "%Ey/%b/%0d%t%0I:%0M:%0S3 %p %Z")
  expect_error(
    easy_format(second %>%
                  decimal %>%
                  digits(1) ),
    "\n\n|base   |decimal | d\u2018igits|\n|:------|:-------|------:|\n|second |TRUE    |      1|\n\n has no corresponding code\n")
})


expect_equal(
  parse_binding_errors("\u2018dog\u2019 \n \u2018cat\u2019"),
  "utils::globalVariables(c('dog', 'cat'))" )

expect_equal(
  format_for_output(easy_format("a", "b") ),
  "ab\n")

expect_equal(
  name_value_to_list(data.frame(name = "a", value = 1)),
  c(a = 1) )

expect_equal(
  modify_option("hello"),
  structure(list(type = "raw", final_code = "hello", variable = NA), .Names = c("type",
                                                                                "final_code", "variable"), class = c("tbl_df", "tbl", "data.frame"
                                                                                ), row.names = c(NA, -1L)) )

expect_equal(
  modify_flag("hello"),
  structure(list(type = "raw", final_code = "hello", variable = TRUE), .Names = c("type",
                                                                                  "final_code", "variable"), class = c("tbl_df", "tbl", "data.frame"
                                                                                  ), row.names = c(NA, -1L)) )

expect_equal(
  function_modify(function() list(a = 1),
                  . %>% stringi::stri_replace_all_fixed("a", "b") ),
  function ()
    list(b = 1) )

# looks like you can't join dataframes with list type variables... had to separate columns
expect_equal(
  modify_function_frame(c("a", "b"), modify_flag)[[1]],
  structure(list(name = c("a", "b"), value = list(function (args,
                                                            value = TRUE)
    args %>% string_table %>% dplyr::mutate(a = value), function (args,
                                                                  value = TRUE)
      args %>% string_table %>% dplyr::mutate(b = value))), .Names = c("name",
                                                                       "value"), class = c("rowwise_df", "tbl_df", "tbl", "data.frame"
                                                                       ), row.names = c(NA, -2L)) [[1]])

expect_equal(
  modify_function_frame(c("a", "b"), modify_flag)[[2]],
  structure(list(name = c("a", "b"), value = list(function (args,
                                                            value = TRUE)
    args %>% string_table %>% dplyr::mutate(a = value), function (args,
                                                                  value = TRUE)
      args %>% string_table %>% dplyr::mutate(b = value))), .Names = c("name",
                                                                       "value"), class = c("rowwise_df", "tbl_df", "tbl", "data.frame"
                                                                       ), row.names = c(NA, -2L)) [[2]])
