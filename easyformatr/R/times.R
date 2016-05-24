#' @include tools.R
NULL

time_info.raw =
  dplyr::bind_rows( dplyr::data_frame(code = 'p', base = 'am_pm'),
                    dplyr::data_frame(code = 'EC', base = 'century', religious = TRUE),
                    dplyr::data_frame(code = 'C', base = 'century'),
                    dplyr::data_frame(code = 'D', base = 'date', month_first = TRUE),
                    dplyr::data_frame(code = 'Ex', base = 'date', local_output = TRUE, religious = TRUE),
                    dplyr::data_frame(code = 'x', base = 'date', local_output = TRUE),
                    dplyr::data_frame(code = 'F', base = 'date'),
                    dplyr::data_frame(code = 'Ec', base = 'datetime', religious = TRUE),
                    dplyr::data_frame(code = 'c', base = 'datetime'),
                    dplyr::data_frame(code = '+', base = 'datetime', with_timezone = TRUE),
                    dplyr::data_frame(code = '0d', base = 'day', roman = TRUE),
                    dplyr::data_frame(code = 'd', base = 'day'),
                    dplyr::data_frame(code = 'e', base = 'day', strip_zeros = TRUE),
                    dplyr::data_frame(code = 'a', base = 'day_of_week', name = TRUE, short = TRUE),
                    dplyr::data_frame(code = 'A', base = 'day_of_week', name = TRUE),
                    dplyr::data_frame(code = 'u', base = 'day_of_week', strip_zeros = TRUE),
                    dplyr::data_frame(code = '0w', base = 'day_of_week', roman = TRUE),
                    dplyr::data_frame(code = 'w', base = 'day_of_week'),
                    dplyr::data_frame(code = 'j', base = 'day_of_year'),
                    dplyr::data_frame(code = '0H', base = 'hour', roman = TRUE),
                    dplyr::data_frame(code = 'H', base = 'hour'),
                    dplyr::data_frame(code = '0I', base = 'hour', twelve = TRUE, roman = TRUE),
                    dplyr::data_frame(code = 'I', base = 'hour', twelve = TRUE),
                    dplyr::data_frame(code = '0M', base = 'minute', roman = TRUE),
                    dplyr::data_frame(code = 'M', base = 'minute'),
                    dplyr::data_frame(code = 'b', base = 'month', name = TRUE, short = TRUE),
                    dplyr::data_frame(code = 'B', base = 'month', name = TRUE),
                    dplyr::data_frame(code = '0m', base = 'month', roman = TRUE),
                    dplyr::data_frame(code = 'm', base = 'month'),
                    dplyr::data_frame(code = 'n', base = 'newline'),
                    dplyr::data_frame(code = '0S', base = 'second', decimal = TRUE),
                    dplyr::data_frame(code = 'S', base = 'second'),
                    dplyr::data_frame(code = 'r', base = 'time', twelve = TRUE),
                    dplyr::data_frame(code = 'R', base = 'time'),
                    dplyr::data_frame(code = 'T', base = 'time', with_seconds = TRUE),
                    dplyr::data_frame(code = 'EX', base = 'time', local_output = TRUE, religious = TRUE, with_seconds = TRUE),
                    dplyr::data_frame(code = 'X', base = 'time', local_output = TRUE, with_seconds = TRUE),
                    dplyr::data_frame(code = 'k', base = 'time', blanks = TRUE),
                    dplyr::data_frame(code = 'l', base = 'time', twelve = TRUE, blanks = TRUE),
                    dplyr::data_frame(code = 's', base = 'timestamp'),
                    dplyr::data_frame(code = 'Z', base = 'timezone'),
                    dplyr::data_frame(code = 'z', base = 'timezone', number = TRUE),
                    dplyr::data_frame(code = '0U', base = 'week_of_year', US = TRUE, roman = TRUE),
                    dplyr::data_frame(code = 'U', base = 'week_of_year', US = TRUE),
                    dplyr::data_frame(code = '0V', base = 'week_of_year', roman = TRUE),
                    dplyr::data_frame(code = 'V', base = 'week_of_year'),
                    dplyr::data_frame(code = '0W', base = 'week_of_year', UK = TRUE, roman = TRUE),
                    dplyr::data_frame(code = 'W', base = 'week_of_year', UK = TRUE),
                    dplyr::data_frame(code = 't', base = 'tab'),
                    dplyr::data_frame(code = '0y', base = 'year', short = TRUE, roman = TRUE),
                    dplyr::data_frame(code = 'y', base = 'year', short = TRUE),
                    dplyr::data_frame(code = 'Y', base = 'year'),
                    dplyr::data_frame(code = 'Ey', base = 'year', short = TRUE, religious = TRUE),
                    dplyr::data_frame(code = 'EY', base = 'year', religious = TRUE),
                    dplyr::data_frame(code = '0S0', base = 'second', digits = 0),
                    dplyr::data_frame(code = '0S1', base = 'second', digits = 1),
                    dplyr::data_frame(code = '0S2', base = 'second', digits = 2),
                    dplyr::data_frame(code = '0S3', base = 'second', digits = 3),
                    dplyr::data_frame(code = '0S4', base = 'second', digits = 4),
                    dplyr::data_frame(code = '0S5', base = 'second', digits = 5),
                    dplyr::data_frame(code = '0S6', base = 'second', digits = 6) )

mutant =
  time_info.raw %>%
  dplyr::select(-base) %>%
  tidyr::gather(variable, value,
                -code) %>%
  dplyr::select(-variable) %>%
  dplyr::filter(value %>% is.na %>% magrittr::not()) %>%
  dplyr::select(-value) %>%
  dplyr::distinct() %>%
  dplyr::mutate(component.mutant = "mutant")

time_info =
  time_info.raw %>%
  dplyr::left_join(mutant) %>%
  dplyr::mutate(component =
                  component.mutant %>%
                  is.na %>%
                  ifelse("base", component.mutant)) %>%
  dplyr::select(-component.mutant)

time_base_object =
  time_info %>%
  dplyr::select(-code) %>%
  dplyr::filter(component == "base") %>%
  dplyr::select(-component) %>%
  dplyr::mutate(type = "time") %>%
  dplyr::group_by(base) %>%
  dplyr::do(value = dplyr::as_data_frame(.)) %>%
  dplyr::rename(name = base)

time_option_object =
  time_info %>%
  Filter(. %>% is.numeric, .) %>%
  names %>%
  modify_function_frame(modify_option)

time_flag_object =
  time_info %>%
  Filter(. %>% is.logical, .) %>%
  names %>%
  modify_function_frame(modify_flag)

time_table_to_final_table = function(time_table)
  if (nrow(time_table) == 0) dplyr::data_frame() else {

    # wipe clean non-relevant fields
    wipe_irrelevant = function(one_base_df) {
      relevant_columns =
        one_base_df %>%
        dplyr::select(base) %>%
        dplyr::slice(1) %>%
        dplyr::left_join(time_info) %>%
        dplyr::select(-code) %>%
        remove_NA_columns() %>%
        names %>%
        c("ID") %>%
        dplyr::intersect(names(one_base_df))

      one_base_df %>%
        dplyr::select(one_of(relevant_columns))
    }

    almost_there =
      time_table %>%
      dplyr::group_by(base) %>%
      dplyr::do( wipe_irrelevant(.) ) %>%
      dplyr::left_join(time_info)

    if (NA %in% almost_there$code)
      almost_there %>%
      dplyr::select(-ID) %>%
      dplyr::filter(code %>% is.na) %>%
      remove_NA_columns %>%
      knitr::kable() %>%
      paste(collapse = "\n") %>%
      paste0("\n\n", . , "\n\n has no corresponding code") %>%
      stop

    almost_there %>%
      dplyr::mutate(final_code = paste0("%", code))
  }
