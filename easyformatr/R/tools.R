utils::globalVariables(c('type', 'ID', 'final_code', '.', 'name', 'code', 'value', 'use_input', 'input_dollar', 'whole_flag', 'before_decimal', 'after_decimal', 'n', 'variable', '.id', '.result', 'base', 'one_of'))

#' @importFrom magrittr %>%
#' @export
#' @title Gives you a command to handle "no visible binding" errors
#' @param string The pasted text of no visible binding errors
#' @description See title
parse_binding_errors = function(string)
  string %>%
  stringi::stri_extract_all_regex("\u2018.*\u2019") %>%
  dplyr::first() %>%
  unique %>%
  stringi::stri_sub(2, -2) %>%
  sprintf("'%s'", .) %>%
  paste(collapse = ", ") %>%
  paste0("utils::globalVariables(c(", . , "))")

#' @title Will paste, omitting all NA values
#' @param ... See paste documentation
#' @param sep See paste documentation
#' @param collapse See paste documentation
#' @description See title
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

#' @title Will paste. Anything pasted to NA will return NA
#' @param ... See paste documentation
#' @param sep See paste documentation
#' @param collapse See paste documentation
#' @description See title
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
#' @title Build objects and dataframes together tidily
#' @param x A list, dataframe, string, or object that can be converted to a string.
#' @description This will recursively bind dataframes together rowwise.
#' Anything that is not a list or a dataframe will be converted a string, doubling all percentage signs,
#' and then to dplyr::0data_frame(base = "raw", raw_value = .)
string_table = function(x) UseMethod("string_table")
#' @export
string_table.default = function(x)
  x %>%
  as.character %>%
  stringi::stri_replace_all_fixed("%", "%%") %>%
  dplyr::data_frame(type = "raw", final_code = .)
#' @export
string_table.easy_format = function(x)
  x %>%
  as.character %>%
  dplyr::data_frame(type = "raw", final_code = .)
#' @export
string_table.data.frame = function(x) x
#' @export
string_table.list = function(x)
  x %>%
  lapply(string_table) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(ID = 1:n())

format_for_output = function(x)
  x %>%
  strwrap %>%
  paste(collapse = "\n") %>%
  paste0("\n")

#' @export
print.easy_format = function(x, ...)
  cat(x %>% format_for_output)

modify_option = function(args, value = NA)
  args %>%
  string_table %>%
  dplyr::mutate(variable = value)

modify_flag = function(args, value = TRUE)
  args %>%
  string_table %>%
  dplyr::mutate(variable = value)

function_modify = function(modify_function, string_function)
  modify_function %>%
  deparse %>%
  paste(collapse = "\n") %>%
  string_function %>%
  parse(text = .) %>%
  eval

modify_function_frame = function(function_names, function_to_modify)
  dplyr::data_frame(name = function_names) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(value =
                  function_to_modify %>%
                  function_modify(. %>% stringi::stri_replace_all_fixed("variable", name)) %>%
                  list)

#' @export
#' @title Convert a dataframe that has a name and a value column to a named list.
#' @param x A dataframe
#' @description Convenience function
name_value_to_list = function(x) x$value %>% stats::setNames(x$name)

#' @export
#' @title Remove columns that are all NA
#' @param x A dataframe or list of vectors.
#' @description Convenience function
#' @examples
#' remove_NA_columns(data.frame(a = NA, b = 1))
remove_NA_columns = function(x) x %>% Filter(. %>% is.na %>% all %>% magrittr::not(), .)

one_value = function(vector)
  (vector %in%
     vector[1]) %>%
  all

#' @export
#' @title Remove columns that are all identical
#' @param x A dataframe or list of vectors.
#' @description Convenience function
#' @examples
#' remove_single_value_columns(data.frame(a = c(1, 1), b = c(1, 2)))
remove_single_value_columns = function(x) x %>% Filter(. %>% one_value %>% magrittr::not(), .)
