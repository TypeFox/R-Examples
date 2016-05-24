#' @include tools.R
NULL

number_base_info =
  dplyr::data_frame(integer = "i",
                    octal = "o",
                    hex = "x",
                    double = "f",
                    scientific = "e",
                    auto = "g",
                    binary = "a") %>%
  tidyr::gather(base, code)

number_flag_info.wide =
  dplyr::data_frame(left_justify = "-",
                    always_sign = "+",
                    prefix_space = " ",
                    zero_pad = "0",
                    hex_prefix = "#",
                    always_decimal = "#",
                    remove_zeros = "#")

number_flag_info =
  number_flag_info.wide %>%
  tidyr::gather(flag, code)

number_option_info.wide =
  dplyr::data_frame(before_decimal = "number of digits before the decimal place",
                    after_decimal = "number of digits after the decimal place",
                    use_input = "sprintf input argument number")

number_option_info =
  number_option_info.wide %>%
  tidyr::gather(option, description)

number_base_object =
  number_base_info %>%
  merge(number_option_info.wide %>% magrittr::inset( , , NA) ) %>%
  merge(number_flag_info.wide %>% magrittr::inset( , , FALSE) ) %>%
  dplyr::mutate(type = "number") %>%
  dplyr::group_by(base) %>%
  dplyr::do(. , value = .) %>%
  dplyr::rename(name = base)

number_flag_object =
  number_flag_info %>%
  magrittr::use_series(flag) %>%
  modify_function_frame(modify_flag)

number_option_object =
  number_option_info %>%
  magrittr::use_series(option) %>%
  modify_function_frame(modify_option)

number_table_to_final_table = function(number_table)

  if (nrow(number_table) == 0) dplyr::data_frame() else {

    flag_names =
      names(number_table) %>%
      dplyr::intersect(number_flag_info$flag)

    flags = if (length(flag_names) == 0) dplyr::data_frame(type = "number") else {
      number_table %>%
        dplyr::select(-code) %>%
        tidyr::gather_("flag", "value", flag_names ) %>%
        dplyr::filter(value) %>%
        dplyr::inner_join(number_flag_info) %>%
        dplyr::group_by(ID) %>%
        dplyr::summarize(whole_flag = code %>% paste(collapse = "") )
    }

    number_table %>%
      dplyr::left_join(flags) %>%
      dplyr::left_join(number_base_info) %>%
      dplyr::mutate(input_dollar = paste_NA_poison(use_input, "$"),
                    final_code = paste_NA_omit("%",
                                               input_dollar,
                                               whole_flag,
                                               before_decimal, ".", after_decimal,
                                               code) ) }
