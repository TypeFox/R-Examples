#' @importFrom magrittr %>%

"string_table.default: no visible binding for global variable '.'
string_table.list: no visible global function definition for 'n'
paste_NA_omit: no visible binding for global variable '.'
paste_NA_omit: no visible global function definition for 'n'
paste_NA_omit: no visible binding for global variable 'variable'
paste_NA_omit: no visible binding for global variable 'value'
paste_NA_omit: no visible binding for global variable '.id'
paste_NA_poison: no visible binding for global variable '.'
paste_NA_poison: no visible binding for global variable '.result'
string_format_: no visible binding for global variable 'value'
string_format_: no visible binding for global variable 'name'
string_format_: no visible binding for global variable 'flag'
string_format_: no visible binding for global variable 'ID'
string_format_: no visible binding for global variable 'base'
string_format_: no visible binding for global variable 'base_code'
string_format_: no visible binding for global variable 'use_input'
string_format_: no visible binding for global variable 'input_dollar'
string_format_: no visible binding for global variable 'whole_flag'
string_format_: no visible binding for global variable 'before_decimal'
string_format_: no visible binding for global variable 'after_decimal'
string_format_: no visible binding for global variable 'raw_value'
string_format_: no visible binding for global variable 'final_format'
string_format_: no visible binding for global variable '.'" %>%
  stringi::stri_extract_all_regex("'.*'") %>%
  dplyr::first() %>%
  stringi::stri_sub(2, -2) %>%
  unique %>%
  utils::globalVariables()

paste_NA_omit = function(... , sep = "", collapse = NULL) {
  data =
    list(...) %>%
    do.call("data_frame" %>% utils::getFromNamespace("dplyr"), . ) %>%
    stats::setNames(1:ncol(.))

  if (".id" %in% names(data))
    stop(".id not allowed as a column name")

  vector =
    data %>%
    dplyr::mutate(.id = 1:n()) %>%
    tidyr::gather(variable, value, -.id) %>%
    dplyr::filter(value %>% is.na %>% `!`) %>%
    dplyr::group_by(.id) %>%
    dplyr::summarize(value = paste(value, collapse = sep)) %>%
    magrittr::use_series(value)

  if (collapse %>% is.null) vector else vector %>% paste(collapse = collapse)
}

paste_NA_poison = function(..., sep = "", collapse = NULL) {

  data =
    list(...) %>%
    do.call("data_frame" %>% utils::getFromNamespace("dplyr"), . )

  if (".result" %in% names(data))
    stop(".result not allowed as a column name")

  vector =
    data %>%
    dplyr::mutate(.result =
                    as.list(.) %>%
                    c(sep = sep) %>%
                    do.call("paste", .),
                  .result =
                    stats::complete.cases(.) %>%
                    ifelse(.result, NA)) %>%
    magrittr::use_series(.result)

  if (collapse %>% is.null) vector else
    if (NA %in% vector) NA else
      vector %>% paste(collapse = collapse)
}

#' @export
#' @title Build string format data into handy dataframe form.
#' @param x A list, dataframe, string, or object that can be converted to a string.
#' @description This is an internal function used for formatting. However, its available
#' to users because tidy data is best data. It will recursively bind dataframes together rowwise,
#' converting any non-list to a string, doubling all percentage signs,
#' and then to dplyr::data_frame(base = "raw", raw_value = .)
#' @examples
#' string_table(list(data.frame(base = "raw", raw_value = "1%%"),
#'                  "1%") )
string_table = function(x) UseMethod("string_table")
#' @export
string_table.default = function(x)
  x %>%
  as.character %>%
  stringi::stri_replace_all_fixed("%", "%%") %>%
  dplyr::data_frame(base = "raw", raw_value = .)
#' @export
string_table.data.frame = function(x) x
#' @export
string_table.list = function(x)
  x %>%
  lapply(string_table) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(ID = 1:n())

#' @export
#' @title Get a list of base string components for building string formats
#' @description Will return information about the base string components that can be used in string_formats.
#' @examples
#' string_base
#' string_format(integer, " ", double)
string_base =
  dplyr::data_frame(integer = "i",
                    octal = "o",
                    hex = "x",
                    double = "f",
                    scientific = "e",
                    auto = "g",
                    binary = "a",
                    string = "s",
                    percent = "%")

#' @export
#' @title Get a list of optional string format flags
#' @description Each flag is associated with a function. This function will append the flag to an input string format.
#'     Flag functions take two values: a format to modify, and the value to modify the flag to
#'     (by default, TRUE). Flags can be looked up by code in the sprintf documentation.
#' @examples
#'     string_format(zero_pad(double))
#'     string_format(zero_pad(double, FALSE))
string_flag =
  dplyr::data_frame(left_justify = "-",
                    always_sign = "+",
                    prefix_space = " ",
                    zero_pad = "0",
                    hex_prefix = "#",
                    always_decimal = "#",
                    remove_zeros = "#")

#' @export
#' @title Get a list of string options
#' @description Each string option is associated with a function. This function takes two arguments:
#'     a string format to be modified and value to modify it with. More information is about each option
#'     is included in the returned table.
#' @examples
#'     string_format(after_decimal(double, 3))
string_option =
  dplyr::data_frame(before_decimal = "number of digits before the decimal place",
                    after_decimal = "number of digits after the decimal place",
                    use_input = "sprintf input argument number")

default_string_base =
  string_base %>%
  tidyr::gather(base, base_code) %>%
  merge(string_option %>% magrittr::inset( , , NA) ) %>%
  merge(string_flag %>% magrittr::inset( , , FALSE) ) %>%
  dplyr::select(-base_code) %>%
  dplyr::group_by(base) %>%
  dplyr::do(string_frame = dplyr::slice(., 1) )

string_mutate_function =
  dplyr::data_frame(variable_name = c(names(string_flag),
                                      names(string_option))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(FUN =
                  "function(args, value = TRUE)
                args %>%
                string_table %>%
                dplyr::mutate(variable_name = value)" %>%
                  stringi::stri_replace_all_fixed("variable_name", variable_name) %>%
                  parse(text = .) %>%
                  eval %>%
                  list)

variable_name = "use_input"
string_list =
  c(with(string_mutate_function,
         FUN %>%
           stats::setNames(variable_name) ),
    with(default_string_base,
         string_frame %>%
           stats::setNames(base) ) )

#' @export
#' @title Build a string format.
#' @description Components will be concatenated separated by sep (lists and vectors will be unlisted).
#'     Raw % signs will be doubled.
#'     Special elements can be included: string bases, flags, and options.
#'     See the documentation in string_base, string_flag, and string_option.
#' @param ... String format pieces
#' @param sep String format separator
#' @examples
#' string_format(integer, " ", double)
#' string_format(integer, double, sep = " ")
string_format = function(..., sep = "")
  lazyeval::lazy_dots(...) %>%
  string_format_(sep = sep)

#' @export
#' @title Standard evaluation version of string_format
#' @param args A list of unevaluated time format pieces
#' @param sep A separator for pasting pieces together.
#' @description Useful for interactively building time formats.
#' @examples
#' base_list = c("double", "integer")
#' string_format_(list(base_list[1], "' '", base_list[2]) )
string_format_ = function(args, sep = "") suppressWarnings(suppressMessages({

  format_frame =
    args %>%
    # format_frame = lazyeval::lazy(double %>% use_input(1)) %>%
    lazyeval::as.lazy_dots() %>%
    lazyeval::lazy_eval(string_list) %>%
    string_table

  format_frame %>%
    tidyr::gather_("name", "value", names(string_flag)) %>%
    dplyr::filter(value) %>%
    dplyr::inner_join(string_flag %>%
                        tidyr::gather(name, flag)) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarize(whole_flag = flag %>% paste(collapse = "") ) %>%
    dplyr::right_join(format_frame) %>%
    dplyr::left_join(string_base %>%
                       tidyr::gather(base, base_code)) %>%
    dplyr::mutate(input_dollar = paste_NA_poison(use_input, "$"),
                  string_format = paste_NA_omit("%",
                                                input_dollar,
                                                whole_flag,
                                                before_decimal, ".", after_decimal,
                                                base_code),
                  final_format =
                    base %>%
                    `==`("raw") %>%
                    ifelse(raw_value, string_format) ) %>%
    magrittr::use_series(final_format) %>%
    as.list %>%
    magrittr::inset2("sep", sep) %>%
    do.call("paste", .)
}))
