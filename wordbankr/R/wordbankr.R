#' @importFrom magrittr "%>%"
NULL

#' Connect to the Wordbank database
#'
#' @param mode A string indicating connection mode: one of \code{"local"},
#'   or \code{"remote"} (defaults to \code{"remote"})
#' @return A \code{src} object which is connection to the Wordbank database
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' wordbank <- connect_to_wordbank()
#' rm(wordbank)
#' }
connect_to_wordbank <- function(mode = "remote") {
  
  assertthat::assert_that(is.element(mode, c("local", "remote")))
  address <- switch(mode,
                    local = "localhost",
                    remote = "server.wordbank.stanford.edu")
  
  src <- dplyr::src_mysql(host = address, dbname = "wordbank",
                          user = "wordbank", password = "wordbank")
  return(src)
}


#' Connect to an instrument's Wordbank table
#'
#' @param src A connection to the Wordbank database
#' @param language A string of the instrument's language (insensitive to case
#'   and whitespace)
#' @param form A string of the instrument's form (insensitive to case and
#'   whitespace)
#' @return A \code{tbl} object containing the instrument's data
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' src <- connect_to_wordbank()
#' eng_ws <- get_instrument_table(src, "english", "ws")
#' rm(src, eng_ws)
#' }
get_instrument_table <- function(src, language, form) {
  table_name <- paste(unlist(c("instruments",
                               stringr::str_split(tolower(language), " "),
                               stringr::str_split(tolower(form), " "))),
                      collapse = "_")
  instrument_table <- dplyr::tbl(src, table_name)
  return(instrument_table)
}


#' Connect to a single Wordbank common table
#'
#' @param src A connection to the Wordbank database
#' @param name A string indicating the name of a common table
#' @return A \code{tbl} object
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' src <- connect_to_wordbank()
#' instruments <- get_common_table(src, "instrument")
#' rm(src, instruments)
#' }
get_common_table <- function(src, name) {
  common_table <- dplyr::tbl(src, paste("common", name, sep = "_"))
  return(common_table)
}


#' Get the Wordbank instruments
#'
#' @return A data frame
#' @inheritParams connect_to_wordbank
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' instruments <- get_instruments()
#' }
#' @export
get_instruments <- function(mode = "remote") {
  
  src <- connect_to_wordbank(mode = mode)
  
  instruments <- get_common_table(src, name = "instrument") %>%
    dplyr::rename_(instrument_id = "id") %>%
    dplyr::collect()
  
  rm(src)
  gc()
  
  return(instruments)
  
}


filter_query <- function(filter_language = NULL, filter_form = NULL,
                         mode = "remote") {
  if (!is.null(filter_language) | !is.null(filter_form)) {
    instruments <- get_instruments(mode = mode)
    if (!is.null(filter_language)) {
      instruments <- instruments %>%
        dplyr::filter_(.dots = list(~language == filter_language))
    }
    if (!is.null(filter_form)) {
      instruments <- instruments %>%
        dplyr::filter_(.dots = list(~form == filter_form))
    }
    instrument_ids <- instruments$instrument_id
    return(sprintf("WHERE instrument_id IN (%s)",
                   paste(instrument_ids, collapse = ", ")))
  } else {
    return("")
  }
}


#' Get the Wordbank by-administration data
#'
#' @param language An optional string specifying which language's
#'   administrations to retrieve.
#' @param form An optional string specifying which form's administrations to
#'   retrieve.
#' @param filter_age A logical indicating whether to filter the administrations
#'   to ones in the valid age range for their instrument
#' @inheritParams connect_to_wordbank
#' @return A data frame where each row is a CDI administration and each column
#'   is a variable about the administration (\code{data_id}, \code{age},
#'   \code{comprehension}, \code{production}), its instrument (\code{language},
#'   \code{form}), or its child (\code{birth_order}, \code{ethnicity},
#'   \code{sex}, \code{mom_ed}).
#'
#' @examples
#' \dontrun{
#' english_ws_admins <- get_administration_data("English", "WS")
#' all_admins <- get_administration_data()
#' }
#' @export
get_administration_data <- function(language = NULL, form = NULL,
                                    filter_age = TRUE, mode = "remote") {
  
  src <- connect_to_wordbank(mode = mode)
  
  mom_ed <- get_common_table(src, "momed") %>%
    dplyr::collect() %>%
    dplyr::rename_(momed_id = "id", momed_level = "level",
                   momed_order = "order") %>%
    dplyr::arrange_("momed_order") %>%
    dplyr::transmute_(momed_id = ~as.numeric(momed_id),
                      mom_ed = ~factor(momed_level, levels = momed_level))
  
  admin_query <- paste(
    "SELECT data_id, age, comprehension, production, language, form,
    birth_order, ethnicity, sex, momed_id, age_min, age_max
    FROM common_administration
    LEFT JOIN common_instrument
    ON common_administration.instrument_id = common_instrument.id
    LEFT JOIN common_child
    ON common_administration.child_id = common_child.id",
    filter_query(language, form, mode = mode),
    sep = "\n")
  
  admins <- dplyr::tbl(src, dplyr::sql(admin_query)) %>%
    dplyr::collect() %>%
    dplyr::mutate_(data_id = ~as.numeric(data_id)) %>%
    dplyr::left_join(mom_ed) %>%
    dplyr::select_("-momed_id") %>%
    dplyr::mutate_(sex = ~factor(sex, levels = c("F", "M", "O"),
                                 labels = c("Female", "Male", "Other")),
                   ethnicity = ~factor(ethnicity,
                                       levels = c("A", "B", "O", "W", "H"),
                                       labels = c("Asian", "Black", "Other",
                                                  "White", "Hispanic")),
                   birth_order = ~factor(birth_order,
                                         levels = c(1, 2, 3, 4, 5, 6, 7, 8),
                                         labels = c("First", "Second", "Third",
                                                    "Fourth", "Fifth", "Sixth",
                                                    "Seventh", "Eighth")))
  
  rm(src)
  gc()
  
  if (filter_age) admins <- admins %>%
    dplyr::filter_(.dots = list(~age >= age_min, ~age <= age_max))
  
  admins <- admins %>%
    dplyr::select_(.dots = list("-age_min", "-age_max"))
  return(admins)
  
}


strip_item_id <- function(item_id) {
  as.numeric(stringr::str_sub(item_id, 6, stringr::str_length(item_id)))
}


#' Get the Wordbank by-item data
#'
#' @param language An optional string specifying which language's items to
#'   retrieve.
#' @param form An optional string specifying which form's items to retrieve.
#' @inheritParams connect_to_wordbank
#' @return A data frame where each row is a CDI item and each column is a
#'   variable about it (\code{language}, \code{form}, \code{type},
#'   \code{lexical_category}, \code{category}, \code{uni_lemma}, \code{item},
#'   \code{definition}, \code{num_item_id}).
#'
#' @examples
#' \dontrun{
#' english_ws_items <- get_item_data("English", "WS")
#' all_items <- get_item_data()
#' }
#' @export
get_item_data <- function(language = NULL, form = NULL, mode = "remote") {
  
  src <- connect_to_wordbank(mode = mode)
  
  item_query <- paste(
    "SELECT item_id, definition, language, form, type, name AS category,
    lexical_category, lexical_class, uni_lemma, complexity_category
    FROM common_iteminfo
    LEFT JOIN common_instrument
    ON common_iteminfo.instrument_id = common_instrument.id
    LEFT JOIN common_category
    ON common_iteminfo.category_id = common_category.id
    LEFT JOIN common_itemmap
    ON common_iteminfo.map_id = common_itemmap.uni_lemma",
    filter_query(language, form, mode = mode),
    sep = "\n")
  
  items <- dplyr::tbl(src, dplyr::sql(item_query)) %>%
    dplyr::collect() %>%
    dplyr::mutate_(num_item_id = ~strip_item_id(item_id))
  
  rm(src)
  gc()
  
  return(items)
  
}


#' Get the Wordbank administration-by-item data
#'
#' @param instrument_language A string of the instrument's language (insensitive
#'   to case and whitespace)
#' @param instrument_form A string of the instrument's form (insensitive to case
#'   and whitespace)
#' @param items A character vector of column names of \code{instrument_table} of
#'   items to extract. If not supplied, defaults to all the columns of
#'   \code{instrument_table}
#' @param administrations Either a logical indicating whether to include
#'   administration data or a data frame of administration data (from
#'   \code{get_administration_data})
#' @param iteminfo Either a logical indicating whether to include item data or a
#'   data frame of item data (from \code{get_item_data})
#' @inheritParams connect_to_wordbank
#' @return A data frame where each row is the result (\code{value}) of a given
#'   item (\code{num_item_id}) for a given administration (\code{data_id})
#'
#' @examples
#' \dontrun{
#' eng_ws_data <- get_instrument_data(instrument_language = "English",
#'                                    instrument_form = "WS",
#'                                    items = c("item_1", "item_42"))
#' }
#' @export
get_instrument_data <- function(instrument_language, instrument_form,
                                items = NULL, administrations = FALSE,
                                iteminfo = FALSE, mode = "remote") {
  
  src <- connect_to_wordbank(mode = mode)
  instrument_table <- get_instrument_table(src, instrument_language,
                                           instrument_form)
  
  if (is.null(items)) {
    columns <- instrument_table$select
    items <- as.character(columns)[2:length(columns)]
  } else {
    assertthat::assert_that(all(items %in% instrument_table$select))
    names(items) <- NULL
  }
  
  if (class(administrations) == "logical" && administrations) {
    administrations <- get_administration_data(instrument_language,
                                               instrument_form,
                                               mode = mode)
  }
  
  if (class(iteminfo) == "logical" && iteminfo) {
    iteminfo <- get_item_data(instrument_language, instrument_form,
                              mode = mode) %>%
      dplyr::select_(.dots = list("-language", "-form"))
  }
  
  instrument_data <- instrument_table %>%
    dplyr::select_(.dots = as.list(c("basetable_ptr_id", items))) %>%
    dplyr::collect() %>%
    dplyr::mutate_(data_id = ~as.numeric(basetable_ptr_id)) %>%
    dplyr::select_("-basetable_ptr_id") %>%
    tidyr::gather_("item_id", "value", items) %>%
    dplyr::mutate_(num_item_id = ~strip_item_id(item_id)) %>%
    dplyr::select_("-item_id")
  
  if ("data.frame" %in% class(administrations)) {
    instrument_data <- dplyr::right_join(instrument_data, administrations,
                                         by = "data_id")
  }
  
  if ("data.frame" %in% class(iteminfo)) {
    instrument_data <- dplyr::right_join(instrument_data, iteminfo,
                                         by = "num_item_id")
  }
  
  rm(src, instrument_table)
  gc()
  
  return(instrument_data)
  
}
