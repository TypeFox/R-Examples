#' @include tools.R times.R numbers.R
NULL

time_object =
  list(base = time_base_object,
       flag = time_flag_object,
       option = time_option_object) %>%
  dplyr::bind_rows(.id = "component")

number_object =
  list(base = number_base_object,
       flag = number_flag_object,
       option = number_option_object) %>%
  dplyr::bind_rows(.id = "component")

object =
  list(time = time_object,
       number = number_object) %>%
  dplyr::bind_rows(.id = "type") %>%
  name_value_to_list

number_info =
  list(base = number_base_info,
       flag = number_flag_info,
       option = number_option_info) %>%
  dplyr::bind_rows(.id = "component")

#' @title Format info
#' @description A huge table of documentation. Best queried through filter_info
#' @export
format_info =
  list(time = time_info,
       number = number_info) %>%
  dplyr::bind_rows(.id = "type")

#' @title Query documentation
#' @description Will return a cleaned markdown table
#' @param filter_statement A logical filter
#' @export
#' @examples
#' filter_info(type == "time" & component == "base")
#' filter_info(type == "number" & component == "base")
#' filter_info(type == "time" & component == "mutant" & base == "second")
#' filter_info(type == "number" & component == "flag")
#' filter_info(type == "number" & component == "option")
filter_info = function(filter_statement)
  lazyeval::lazy(filter_statement) %>%
  filter_info_

#' @title Standard evaluation version of filter_info
#' @description Will return a cleaned markdown table
#' @param filter_expression An unevaluated logical filter
#' @export
#' @examples
#' filter_info_("type == 'time' & component == 'base'")
filter_info_ = function(filter_expression) {
  result =
    format_info %>%
    dplyr::filter_(filter_expression) %>%
    remove_NA_columns

  if (nrow(result) == 0)
    stop("No information matches those conditions") else
      if (nrow(result) == 1)
        result %>%
    knitr::kable() else
      result %>%
    remove_single_value_columns %>%
    knitr::kable()
}

#' @title Standard evaluation version of easy_format
#' @description See title
#' @param args An unevaluated expression, or list of expressions
#' @param sep A separator, defaults to ""
#' @export
#' @examples
#' easy_format_(c("hour", "minute", "second"))
easy_format_ = function(args, sep = "") suppressMessages({
  format_table =
    args %>%
    lazyeval::as.lazy_dots() %>%
    lazyeval::lazy_eval(object) %>%
    string_table()

  dplyr::bind_rows(format_table %>%
                     dplyr::filter(type == "time") %>%
                     time_table_to_final_table,
                   format_table %>%
                     dplyr::filter(type == "number") %>%
                     number_table_to_final_table,
                   format_table %>%
                     dplyr::filter(type == "raw") ) %>%
    dplyr::select(ID, final_code) %>%
    dplyr::arrange(ID) %>%
    magrittr::use_series(final_code) %>%
    paste(collapse = sep) %>%
    structure(class = "easy_format")
})

#' @title Easily build format strings
#' @description Use `filter_info` to find more information about format components
#' @param ... A mixture of time or number components, including bases and modifications, and strings.
#' @param sep A separator, defaults to ""
#' @export
#' @examples
#' easy_format(year, month, day, integer, octal, double)
#' easy_format(decimal(second) )
#' easy_format(before_decimal(double, 3) )
#' easy_format(month,
#'             roman(list(day,
#'                        minute) ) )
easy_format = function(..., sep = "")
  lazyeval::lazy_dots(...) %>%
    easy_format_(sep = sep)

