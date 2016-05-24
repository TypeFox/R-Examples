utils::globalVariables(c(".", "base", "code", "variable", "value"))

#' @importFrom magrittr %>%
#' @export
#' @title Build a time format expression
#' @description See vignette
#' @param ... Time format pieces
time_format = function(...)
  lazyeval::lazy_dots(...) %>%
  time_format_

time_code = structure(list(code = c("p", "EC", "C", "D", "Ex", "x", "F",
                                      "Ec", "c", "+", "0d", "d", "e", "a", "A", "u", "0w", "w", "j",
                                      "0H", "H", "0I", "I", "0M", "M", "b", "B", "0m", "m", "n", "%",
                                      "0S0", "0S1", "0S2", "0S3", "0S4", "0S5", "0S6", "0S", "S", "r",
                                      "R", "T", "EX", "X", "k", "l", "s", "Z", "z", "0U", "U", "0V",
                                      "V", "0W", "W", "t", "0y", "y", "Y", "Ey", "EY"),
                           base = c("am_pm",
                                    "century", "century", "date", "date", "date", "date", "datetime",
                                    "datetime", "datetime", "day", "day", "day", "day_of_week", "day_of_week",
                                    "day_of_week", "day_of_week", "day_of_week", "day_of_year", "hour",
                                    "hour", "hour", "hour", "minute", "minute", "month", "month",
                                    "month", "month", "newline", "percent", "second", "second", "second",
                                    "second", "second", "second", "second", "second", "second", "time",
                                    "time", "time", "time", "time", "time", "time", "timestamp",
                                    "timezone", "timezone", "week_of_year", "week_of_year", "week_of_year",
                                    "week_of_year", "week_of_year", "week_of_year", "tab",
                                    "year", "year", "year", "year", "year"),
                           name = c(NA, NA, NA,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE, NA, NA, NA,
                                    NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE, NA, NA, NA, NA, NA, NA,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           short = c(NA,
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA,
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, NA, NA, NA,
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE, NA, TRUE,
                                     NA),
                           strip_zeros = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, TRUE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA),
                           twelve = c(NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, TRUE, NA, NA, NA, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA),
                           US = c(NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           UK = c(NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                  NA, TRUE, TRUE, NA, NA, NA, NA, NA, NA),
                           blanks = c(NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           roman = c(NA,
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, NA, NA, NA,
                                     TRUE, NA, NA, TRUE, NA, TRUE, NA, TRUE, NA, NA, NA, TRUE, NA,
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                     NA, NA, NA, NA, NA, TRUE, NA, TRUE, NA, TRUE, NA, NA, TRUE, NA,
                                     NA, NA, NA),
                           local_output = c(NA, NA, NA, NA, TRUE, TRUE, NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                            NA, NA, NA, NA, TRUE, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA),
                           religious = c(NA, TRUE, NA,
                                         NA, TRUE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                         NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, NA, NA, NA,
                                         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE),
                           decimal = c(NA,
                                       NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                       NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                       NA, NA, NA, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                       NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           digits = c(NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0L, 1L,
                                      2L, 3L, 4L, 5L, 6L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           month_first = c(NA,
                                           NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           with_timezone = c(NA,
                                             NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, NA, NA, NA, NA,
                                             NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                             NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                             NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                           with_seconds = c(NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, TRUE, TRUE, NA, NA,
                                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)),
                      .Names = c("code",
                                 "base", "name", "short", "strip_zeros", "twelve", "US", "UK",
                                 "blanks", "roman", "local_output", "religious", "decimal", "digits",
                                 "month_first", "with_timezone", "with_seconds"),
                      class = "data.frame",
                      row.names = c(NA,
                                    -62L))


base_list =
  time_code %>%
  dplyr::select(-code) %>%
  dplyr::group_by(base) %>%
  dplyr::do(time_frame =
              dplyr::slice(., 1) %>%
              magrittr::inset( , -1, NA))

make_it_a_code =
  "name = function(time_frame, value = TRUE)
time_frame %>%
dplyr::mutate(name = value)"

#' @importFrom stringi stri_replace_all_fixed
make_it_a_list =
  dplyr::data_frame(name = names(time_code %>%
                                   dplyr::select(-code, -base))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(FUN =
                  make_it_a_code %>%
                  stri_replace_all_fixed("name", name) %>%
                  parse(text = .) %>%
                  eval %>%
                  list)

get_string = function(test_frame)
  suppressMessages({
    if (test_frame %>% is.data.frame %>% `!`) {
      test_frame

    } else {
      query =
        test_frame %>%
        dplyr::inner_join(time_code)

      if (nrow(query) > 0) {
        query %>%
          magrittr::use_series(code) %>%
          paste0("%", .)

      } else {
        test_frame %>%
          clean_table %>%
          paste(collapse = "\n") %>%
          paste0("\n", . , "\n\n has no corresponding code") %>%
          stop } } })


time_list =
  c(with(make_it_a_list,
         FUN %>% setNames(name) ),
    with(base_list,
         time_frame %>% setNames(base)))


paste_00 = function(...) paste(sep = "", collapse = "", ...)

#' @export
#' @title Standard evaluation version of time_format
#' @param args A list of unevaluated time format pieces
#' @description See vignette for lazyeval
time_format_ = function(args)
  args %>%
  lazyeval::as.lazy_dots() %>%
  lazyeval::lazy_eval(time_list) %>%
  lapply(get_string) %>%
  paste_00

#' @export
#' @title Get a list of base units of time for building time formats
#' @description See vignette
time_bases = function()
  suppressMessages({
    time_code %>%
      tidyr::gather(variable, value,
                    -code, -base) %>%
      dplyr::filter(value %>% is.na %>% `!`) %>%
      dplyr::anti_join(time_code, .) %>%
      dplyr::select(code, base) %>%
      dplyr::arrange(base) %>%
      knitr::kable() })

#' @export
#' @title Get a list of optional modifications to a unit of time
#' @description See vignette
#' @param specific_base A base unit of time to explore (see vignette)
base_options = function(specific_base)
  time_code %>%
  dplyr::filter(base == specific_base) %>%
  clean_table

clean_table = function(table) {
  answer =
    table %>%
    Filter(. %>% is.na %>% all %>% `!`, . )
  answer[is.na(answer)] = ""
  answer[answer == TRUE] = "1"
  answer %>% knitr::kable()
}
